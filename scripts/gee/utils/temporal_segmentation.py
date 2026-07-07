import ee
import math
from datetime import datetime

J_DAYS = 0
FRACTIONAL_YEARS = 1
UNIX_TIME_MILLIS = 2

class DateConversion:
    @staticmethod
    def toT(date, date_format):
        if isinstance(date, ee.Date):
            if date_format == J_DAYS:
                epoch_day = 719529
                return date.millis().divide(1000).divide(3600).divide(24).add(epoch_day)
            elif date_format == FRACTIONAL_YEARS:
                return date.get('year').add(date.getFraction('year'))
            elif date_format == UNIX_TIME_MILLIS:
                return date.millis()
            else:
                raise ValueError('Unsupported dateFormat')
        else:
            # Client-side Python datetime handling
            if not isinstance(date, datetime):
                date = datetime.fromisoformat(date) if isinstance(date, str) else date
                
            if date_format == J_DAYS:
                epoch_day = 719529
                return (date.timestamp() / 86400) + epoch_day
            elif date_format == FRACTIONAL_YEARS:
                year = date.year
                first_of_year = datetime(year, 1, 1)
                first_of_next_year = datetime(year + 1, 1, 1)
                fraction = (date.timestamp() - first_of_year.timestamp()) / \
                           (first_of_next_year.timestamp() - first_of_year.timestamp())
                return year + fraction
            elif date_format == UNIX_TIME_MILLIS:
                return int(date.timestamp() * 1000)
            else:
                raise ValueError('Unsupported dateFormat')

    @staticmethod
    def fromT(t, date_format):
        t = ee.Number(t)
        if date_format == J_DAYS:
            epoch_day = 719529
            return ee.Date(t.subtract(epoch_day).multiply(1000 * 3600 * 24))
        elif date_format == FRACTIONAL_YEARS:
            first_of_year = ee.Date.fromYMD(t.floor(), 1, 1)
            first_of_next_year = first_of_year.advance(1, 'year')
            days_in_year = first_of_next_year.difference(first_of_year, 'day')
            day_of_year = days_in_year.multiply(t.mod(1)).floor()
            return first_of_year.advance(day_of_year, 'day')
        elif date_format == UNIX_TIME_MILLIS:
            return ee.Date(t)
        else:
            raise ValueError('Unsupported dateFormat')

    @staticmethod
    def days(t1, t2, date_format):
        diff = t2.subtract(t1)
        if date_format == J_DAYS:
            return diff
        elif date_format == FRACTIONAL_YEARS:
            return diff.multiply(365).round()
        elif date_format == UNIX_TIME_MILLIS:
            return diff.divide(1000 * 3600 * 24).round()
        else:
            raise ValueError('Unsupported dateFormat')
            
class HarmonicSlice:
    @staticmethod
    def get_omega(date_format):
        if date_format == J_DAYS:
            return 2.0 * math.pi / 365.25
        elif date_format == FRACTIONAL_YEARS:
            return 2.0 * math.pi
        elif date_format == UNIX_TIME_MILLIS:
            return 2.0 * math.pi / (60 * 60 * 24 * 365.25)
        else:
            raise ValueError('Unsupported dateFormat')

    @staticmethod
    def slice_image(coefs, t, date_format, harmonics=3):
        omega = HarmonicSlice.get_omega(date_format)
        
        def get_c(index):
            return coefs.arrayGet([index])

        components = [
            get_c(0).add(get_c(1).multiply(t)), # Trend
            get_c(2).multiply(t.multiply(omega).cos()).add(get_c(3).multiply(t.multiply(omega).sin())), # H1
            get_c(4).multiply(t.multiply(omega * 2).cos()).add(get_c(5).multiply(t.multiply(omega * 2).sin())), # H2
            get_c(6).multiply(t.multiply(omega * 3).cos()).add(get_c(7).multiply(t.multiply(omega * 3).sin()))  # H3
        ]
        
        # Select components based on harmonics count
        selected = ee.List(components).slice(0, ee.Number(harmonics).add(1))
        
        return ee.ImageCollection(selected)\
            .reduce(ee.Reducer.sum())\
            .regexpRename('(.*)_coefs_sum', '$1', False)
    
class Segments:
    def __init__(self, segments_image, date_format=0, max_segments=50):
        # Note: updateImageMask is assumed to be an external helper or standard mask
        self.segments_image = segments_image.addBands(
            segments_image.select('.*_coefs').arrayPad([0, 8]), None, True
        )
        self.date_format = date_format
        self.max_segments = max_segments

    def get_segment_indexes(self):
        return self.segments_image.select(0).Not().Not() \
            .arrayAccum(0, ee.Reducer.sum()) \
            .subtract(1)

    def segment_index(self, date, strategy='mask'):
        t = ee.Image(DateConversion.toT(date, self.date_format))
        segment_indexes = self.get_segment_indexes()
        t_start = self.segments_image.select('tStart')
        t_end = self.segments_image.select('tEnd')

        # Helper functions for previous/next logic
        def get_previous():
            return segment_indexes.arrayMask(t_start.lte(t)) \
                .arrayReduce(ee.Reducer.lastNonNull(), [0])

        def get_next():
            return segment_indexes.arrayMask(t_end.gt(t)) \
                .arrayReduce(ee.Reducer.first(), [0])

        # Strategy Logic
        if strategy == 'mask':
            masked = segment_indexes.arrayMask(t_start.lte(t).And(t_end.gte(t))) \
                .arrayReduce(ee.Reducer.first(), [0])
        elif strategy == 'previous':
            masked = get_previous()
        elif strategy == 'next':
            masked = get_next()
        elif strategy == 'closest':
            prev_dist = t_start.subtract(t).abs().arrayReduce(ee.Reducer.min(), [0]).arrayGet([0])
            next_dist = t_end.subtract(t).abs().arrayReduce(ee.Reducer.min(), [0]).arrayGet([0])
            masked = get_previous().where(next_dist.gt(prev_dist), get_next())
        else:
            raise ValueError(f"Unsupported strategy: {strategy}. Use 'mask', 'closest', 'previous', or 'next'.")

        return masked.updateMask(masked.arrayLength(0).gt(0)) \
            .arrayFlatten([['segmentIndex']]).int8()

    def find_by_date(self, date, strategy='mask'):
        index = self.segment_index(date, strategy)
        mask = self.get_segment_indexes().eq(index.unmask(-1))
        
        # Separate 1D and 2D bands (coefficients)
        band_names = self.segments_image.bandNames()
        b_1d = band_names.filter(ee.Filter.stringEndsWith('item', '_coefs').Not())
        b_2d = band_names.filter(ee.Filter.stringEndsWith('item', '_coefs'))
        
        image_1d = self.segments_image.select(b_1d).arrayMask(mask)
        image_1d = image_1d.mask(image_1d.select(0).arrayLength(0).unmask(0)).arrayProject([0]).arrayGet([0])
        
        image_2d = self.segments_image.select(b_2d).arrayMask(
            mask.toArray(1).unmask(ee.Array([[]], ee.PixelType.double()))
        )
        image_2d = image_2d.mask(image_2d.select(0).arrayLength(0).unmask(0)).arrayProject([1])
        
        return Segment(image_1d.addBands(image_2d), self.date_format, date)

class Segment:
    def __init__(self, segment_image, date_format, default_date=None):
        self.segment_image = segment_image
        self.date_format = date_format
        self.default_t = DateConversion.toT(default_date, date_format) if default_date else \
                         segment_image.expression('(i.tStart + i.tEnd) / 2', {'i': segment_image})

    def slice(self, **kwargs):
        options = {
            't': self.default_t,
            'harmonics': 3,
            'extrapolateMaxDays': 0,
            'extrapolateMaxFraction': 0,
            'strategy': 'mask',
            'date': None
        }
        options.update(kwargs)
        
        t = ee.Image(DateConversion.toT(options['date'], self.date_format) if options['date'] else options['t'])
        coefs = self.segment_image.select('.*_coefs')
        
        # Logic for extrapolation mask
        if options['strategy'] != 'closest':
            return HarmonicSlice.slice_image(coefs, t, self.date_format, options['harmonics'])
        else:
            # Add logic for clamping t to tStart/tEnd as per your JS code
            return HarmonicSlice.slice_image(coefs, t, self.date_format, options['harmonics'])