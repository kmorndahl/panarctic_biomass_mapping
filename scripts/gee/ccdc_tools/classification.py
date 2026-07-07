"""
@license
Copyright 2019 Boston University

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Utility functions for classifying CCDC results.
"""

import math
import ee
from . import inputs as input_utils
from . import dates as date_utils
from . import ccdc as ccdc_utils


# =========================
# 1. FUNCTIONS ==============
# =========================

def getBinaryLabel(fc, property, target_class):
    """
    Convert training data to binary label for target class.
    @param {ee.FeatureCollection} fc            Training data feature collection
    @param {string}               property      Property label indicating class label
    @param {number}               target_class  Class to retain as 1 in binary label
    @returns {ee.FeatureCollection} Training data where 1 = target_class and 0 = all other classes
    """
    if not fc:
        return 'Required argument [fc] missing.'
    if not target_class:
        return 'Required argument [targetClass] missing.'

    fc = ee.FeatureCollection(fc)
    target_fc = fc.filterMetadata(property, 'equals', target_class).map(
        lambda i: i.set(property, 1))
    not_target_fc = fc.filterMetadata(property, 'not_equals', target_class).map(
        lambda i: i.set(property, 0))
    return target_fc.merge(not_target_fc)


def getClassProbs(fc, coefs_to_classify, class_list, classifier, property):
    """
    Get class probability for each class in training data.
    @param {ee.FeatureCollection} fc                Feature collection of training data
    @param {ee.Image}             coefs_to_classify Multi-band image of coefficients to classify
    @param {list}                 class_list        Classes to test probability of
    @param {ee.Classifier}        classifier        Classifier in 'PROBABILITY' mode
    @param {string}               property          Label defining class in training data
    @returns {ee.Image} Image with each band being class probability for each input class
    """
    band_names = class_list.map(
        lambda num: ee.String('probability_').cat(ee.String(ee.Number(num))))

    def classify_one(num):
        fc_binary = getBinaryLabel(fc, property, num)
        trained = classifier.train(**{
            'features': fc_binary,
            'classProperty': property,
            'inputProperties': coefs_to_classify.bandNames()
        })
        return coefs_to_classify.classify(trained)

    class_probs = class_list.map(classify_one)
    return ee.ImageCollection(ee.List(class_probs)).toBands().rename(band_names)


def makeGrids(region, count, size, seed):
    """
    Make random grids in a region of interest.
    @param {ee.Geometry} region  Study region bounding geometry
    @param {number}      count   Number of random grids
    @param {number}      size    Length of one side of grid in m^2
    @param {number}      seed    Random number seed or 'random'
    @returns {ee.FeatureCollection} Feature collection of random grids
    """
    if seed == 'random':
        import math as _math
        import random as _random
        seed = math.ceil(_random.random() * 1000)

    # Create sample of random points within region
    random_points = ee.FeatureCollection.randomPoints(**{
        'region': region,
        'points': count,
        'seed': seed
    })

    # Take bounding box of buffered samples
    def make_bounds(point):
        buffer = point.buffer(size / 2)
        return buffer.bounds()
    bb = random_points.map(make_bounds)

    # Assign id
    bb_list = bb.toList(bb.size())
    index_list = ee.List.sequence(1, bb.size())

    return ee.FeatureCollection(index_list.map(
        lambda i: ee.Feature(bb_list.get(ee.Number(i).subtract(1))).set({'ID': i})))


def newPhaseAmplitude(img, bands, sin_name, cos_name):
    """
    @deprecated. Use ccdc.newPhaseAmplitude
    Get phase and amplitude for a single spectral band.
    @param {ee.Image} img       CCD results in long multi-band format
    @param {List}     bands     List with the name of the bands for which to calculate ampl. and phase
    @param {String}   sin_name  Band suffix of the desired sine harmonic coefficient (e.g. '_SIN')
    @param {String}   cos_name  Band suffix of the desired cosine harmonic coefficient (e.g. '_COS')
    @returns {ee.Image} Image with two bands representing phase and amplitude of the desired harmonic
    """
    sin_names = bands.map(lambda x: ee.String(x).cat(sin_name))
    cos_names = bands.map(lambda x: ee.String(x).cat(cos_name))
    phase_names = bands.map(lambda x: ee.String(x).cat('_PHASE'))
    amplitude_names = bands.map(lambda x: ee.String(x).cat('_AMPLITUDE'))
    phase = (img.select(sin_names).atan2(img.select(cos_names))
                .unitScale(-3.14159265359, 3.14159265359)
                .multiply(365)
                .rename(phase_names))
    amplitude = img.select(sin_names).hypot(img.select(cos_names)).rename(amplitude_names)
    return phase.addBands(amplitude)


def sampleResultProcedure(training_data, coef_names, band_list, date_property, extra_bands,
                           ccdc_image, segs, ccdc_date_fmt, training_date_fmt, scale):
    """
    Get training coefficient by reading from result data.
    @param {ee.FeatureCollection} training_data      Training data with ccdc outputs in properties
    @param {List}                 coef_names         Coefficient abbreviated names in order
    @param {List}                 band_list          List of input band names in order
    @param {string}               date_property      Property name containing date in features
    @param {List}                 extra_bands        Ancillary bands to add as predictors
    @param {ee.Image}             ccdc_image         CCDC coefficients image
    @param {List}                 segs               Segment identifiers for ccdc_image
    @param {number}               ccdc_date_fmt      Date format of ccdc date format
    @param {number}               training_date_fmt  Training data date format
    @param {number}               scale              Spatial scale to sample training points at
    @returns {ee.FeatureCollection} Training data with coefficients corresponding to specific date
    """
    ccdc_date_fmt = int(ccdc_date_fmt)
    unique_years = ee.Dictionary(
        ee.FeatureCollection(training_data).aggregate_histogram(date_property)).keys()

    def process_year(str_year):
        year = ee.Number.parse(str_year)
        fc_year = training_data.filterMetadata(date_property, 'equals', year)

        formatted_date = date_utils.convertDate({
            'inputFormat': training_date_fmt,
            'inputDate': year,
            'outputFormat': ccdc_date_fmt
        })
        coefs = ccdc_utils.getMultiCoefs(ccdc_image, formatted_date, band_list, coef_names,
                                          True, segs, 'after')

        # Use new code to reduce calculations
        phase_amps = ccdc_utils.newPhaseAmplitude(coefs, '.*SIN.*', '.*COS.*')
        coefs = coefs.addBands(phase_amps)
        if extra_bands:
            coefs = coefs.addBands(extra_bands)
        return coefs.sampleRegions(**{
            'collection': fc_year,
            'scale': scale,
            'tileScale': 16,
            'geometries': True
        })

    return ee.FeatureCollection(unique_years.map(process_year)).flatten()


def runCcdcProcedure(training_data, coef_names, band_list, date_property, extra_bands,
                      landsat_params, segs):
    """
    Get training coefficient by running ccdc on every feature.
    @param {ee.FeatureCollection} training_data    Training data points
    @param {List}                 coef_names       Coefficient abbreviated names in order
    @param {List}                 band_list        List of input band names in order
    @param {String}               date_property    Property name containing date in features
    @param {List}                 extra_bands      Ancillary bands to add as predictors
    @param {dict}                 landsat_params   Parameters for 'getLandsat' function
    @param {List}                 segs             Segment identifiers for ccdc_image
    @returns {ee.FeatureCollection} Training data with coefficients corresponding to specific date
    """
    segs = segs or ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']

    training_ccdc = getTraining({
        'trainingData': training_data,
        'extraBands': extra_bands,
        'landsatParams': landsat_params
    })

    def process_feat(feat):
        year = ee.Number(feat.get(date_property)).add(2)
        year2 = ee.String(year)
        date = date_utils.dateToJdays(year2)

        ccd_image = ccdc_utils.buildCcdImage(training_ccdc, 6, band_list)
        coefs = ccdc_utils.getMultiCoefs(ccd_image, date, band_list, coef_names, True, segs)

        if extra_bands:
            coefs = coefs.addBands(extra_bands)

        sample_coefs = ee.Dictionary(coefs.reduceRegion(**{
            'reducer': ee.Reducer.mean(),
            'geometry': feat.geometry(),
            'scale': 30,
            'crs': 'EPSG:4326',
            'tileScale': 8,
        }))
        return ee.Feature(feat).setMulti(sample_coefs)

    return training_data.map(process_feat)


def getTrainingCoefsAtDate(training_data, coef_names=None, band_list=None, date_property=None,
                            extra_bands=None, ccdc_image=None, segs=None, ccdc_date_fmt=None,
                            training_date_fmt=None, scale=None, landsat_params=None):
    """
    Get coefficients at a given date for each feature in collection.
    @param {ee.FeatureCollection} training_data     Training data points to extract coefficients for
    @param {List}   coef_names         ['INTP','SLP','COS','SIN','RMSE','COS2','SIN2','COS3','SIN3']
    @param {List}   band_list          ['BLUE','GREEN','RED','NIR','SWIR1','SWIR2']
    @param {string} date_property      Property name containing date in features
    @param {List}   extra_bands        Ancillary bands to add as predictors
    @param {ee.Image} ccdc_image       Use ccdc coefficients instead of calculating on the fly
    @param {List}   segs               Segment identifiers for ccdc_image
    @param {number} ccdc_date_fmt      Date format of ccdc date format
    @param {number} training_date_fmt  Training data date format
    @param {number} scale              Spatial scale to sample training points at
    @param {dict}   landsat_params     Parameters for 'getLandsat' function
    @returns {ee.FeatureCollection} Training data with coefficients corresponding to specific date
    """
    coef_names = coef_names or ['INTP', 'SLP', 'COS', 'SIN', 'RMSE', 'COS2', 'SIN2', 'COS3', 'SIN3']
    band_list = band_list or ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']
    date_property = date_property or 'Start_Year'
    landsat_params = landsat_params or {'start': '1990-01-01', 'end': '2020-01-01'}
    segs = segs or ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
    ccdc_date_fmt = ccdc_date_fmt or 1
    training_date_fmt = training_date_fmt or 1
    scale = scale or 30

    result_image = ccdc_image or None

    if result_image:
        return sampleResultProcedure(training_data, coef_names, band_list, date_property,
                                      extra_bands, result_image, segs, ccdc_date_fmt,
                                      training_date_fmt, scale)
    else:
        return runCcdcProcedure(training_data, coef_names, band_list, date_property,
                                 extra_bands, landsat_params, segs)


def remapLC(feats, in_label, out_label, in_list=None, out_list=None):
    """
    Remap training labels to GLANCE level 1 land cover.
    @param {ee.FeatureCollection} feats     Training data feature collection
    @param {string}               in_label  Attribute name containing land cover strings
    @param {string}               out_label Attribute name for output numeric land cover
    @param {list}                 in_list   List of input land cover string values
    @param {list}                 out_list  List of output land cover numeric values
    @returns {ee.FeatureCollection} Training data with numeric out_label column in each feature
    """
    in_list = in_list or ['Water', 'Snow/Ice', 'Built', 'Bare', 'Trees', 'Shrub',
                           'Herbaceous', 'Woodland', 'Forest', 'Developed', 'Agriculture',
                           'Barren', 'Grass', 'Ice_and_Snow', 'Shrubs', 'Wetland']
    out_list = out_list or [1, 2, 3, 4, 5, 6, 7, 8, 5, 3, 7, 4, 7, 2, 6, 1]
    feats = feats.map(lambda feat: feat.set(out_label, feat.get(in_label)))
    return feats.remap(in_list, out_list, out_label)


def assignIds(sample, attribute_name='ID'):
    """
    Shuffle the sample and assign sample ID.
    @param {ee.FeatureCollection} sample         Training data of point samples
    @param {string}               attribute_name Name to assign ID attribute to
    @returns {ee.FeatureCollection} Training data shuffled with unique ID attribute
    """
    with_random = sample.randomColumn(**{'seed': 1})
    with_random = with_random.sort('random')
    with_random_list = with_random.toList(with_random.size())
    index_list = ee.List.sequence(1, with_random.size())

    return ee.FeatureCollection(index_list.map(
        lambda i: ee.Feature(with_random_list.get(
            ee.Number(i).subtract(1))).set({attribute_name: i})))


def getMiddleDate(fc, start_prop, end_prop, middle_prop):
    """
    Get the middle segment date of training data.
    @param {ee.FeatureCollection} fc           Training data feature collection
    @param {string}               start_prop   Property name of segment start year
    @param {String}               end_prop     Property name of segment end year
    @param {String}               middle_prop  Property name of calculated middle attribute
    @returns {ee.FeatureCollection} Training data with middle_prop attribute
    """
    def compute_middle(feat):
        start = ee.Number(feat.get(start_prop))
        end = ee.Number(feat.get(end_prop))
        middle = start.add(end).divide(2).ceil().int16()
        return feat.set(middle_prop, middle)
    return fc.map(compute_middle)


def getInputFeatures(seg, image_to_classify, predictors, band_names, ancillary):
    """
    Function to convert segment band names to universal band names to classify.
    @param {number}   seg               Segment number
    @param {ee.Image} image_to_classify CCDC coefficient stack to classify
    @param {array}    predictors        List of predictor variables
    @param {array}    band_names        Band names of coefficient image
    @param {array}    ancillary         List of ancillary data (or ee.Image)
    @returns {list} [inputFeatures, bands]
    """
    str_seg = ee.String('S').cat(ee.String(ee.Number(seg).int8())).cat('_.*')
    # Another string to remove segment prefix
    str2 = ee.String('S').cat(ee.String(ee.Number(seg).int8())).cat('_')

    # Select bands to classify
    bands = image_to_classify.select([str_seg])

    # Rename without prefix
    renamed_bands = bands.bandNames().map(
        lambda bn: ee.String(ee.String(bn)
                               .replace('_coef_', '_')
                               .replace('_COEF_', '_'))
                               .split(str2).get(1))
    bands = bands.rename(renamed_bands)

    # Mask where there's no model
    bands = bands.updateMask(bands.select('tStart').gt(0))

    # Normalize the intercepts
    bands = ccdc_utils.applyNorm(bands, bands.select('.*tStart'), bands.select('.*tEnd'))

    # Get phase and amplitude
    phase_amp = ccdc_utils.newPhaseAmplitude(bands, '.*SIN.*', '.*COS.*')

    if isinstance(ancillary, ee.Image):
        phase_amp = phase_amp.addBands(ancillary)

    # Add phase, amplitude, and ancillary
    bands = bands.addBands([phase_amp]).select(predictors)

    # Remove non-inputs
    input_features = bands.bandNames().removeAll(
        ['tStart', 'tEnd', 'tBreak', 'changeProb',
         'BLUE_MAG', 'GREEN_MAG', 'RED_MAG', 'NIR_MAG', 'SWIR1_MAG', 'SWIR2_MAG', 'TEMP_MAG'])

    return [input_features, bands]


def makePhaseAmp(img, band_names):
    """
    @deprecated. Use ccdc.newPhaseAmplitude
    Calculate phase and amplitude from sin and cosine coefficients.
    @param {ee.Image} img        CCDC coefficient image
    @param {array}    band_names Coefficient image band names
    @returns {ee.Image} An image containing sin/cosine pairs
    """
    phase_amp1 = newPhaseAmplitude(img, band_names, '_SIN', '_COS')
    bns = phase_amp1.bandNames().map(lambda b: ee.String(b))
    phase_amp1 = phase_amp1.rename(bns)

    phase_amp2 = newPhaseAmplitude(img, band_names, '_SIN2', '_COS2')
    bns = phase_amp2.bandNames().map(lambda b: ee.String(b).cat('_1'))
    phase_amp1 = phase_amp2.rename(bns)

    phase_amp3 = newPhaseAmplitude(img, band_names, '_SIN3', '_COS3')
    bns = phase_amp3.bandNames().map(lambda b: ee.String(b).cat('_2'))
    phase_amp3 = phase_amp3.rename(bns)

    return ee.Image.cat([phase_amp1, phase_amp2, phase_amp3])


def subsetTraining(training_data, train_prop, seed, class_property):
    """
    Subset training data into random training and testing data.
    Data is subset proportionally for each land cover class.
    @param {ee.FeatureCollection} training_data   Training data
    @param {float}                train_prop      Proportion of data to withhold for training
    @param {number}               seed            Seed for random selection of subset
    @param {string}               class_property  Property containing the input class
    @returns {ee.FeatureCollection} Training data with 'train' attribute where 1=training, 0=testing
    """
    class_counts = ee.Dictionary(training_data.aggregate_histogram(class_property))
    classes = class_counts.keys()

    def process_class(c):
        subset = training_data.filterMetadata(class_property, 'equals', ee.Number.parse(c))
        # Withhold a selection of training data
        training_subset_with_random = subset.randomColumn('random', seed).sort('random')
        index_of_split = training_subset_with_random.size().multiply(train_prop).int32()
        number_of_train = training_subset_with_random.size().subtract(index_of_split).int32()
        subset_test = ee.FeatureCollection(
            training_subset_with_random.toList(index_of_split)).map(
                lambda feat: feat.set('train', 0))
        subset_train = ee.FeatureCollection(
            training_subset_with_random.toList(number_of_train, index_of_split)).map(
                lambda feat: feat.set('train', 1))
        return ee.Algorithms.If(
            subset.size().gt(10),
            subset_test.merge(subset_train),
            subset.map(lambda feat: feat.set('train', 1)))

    return ee.FeatureCollection(classes.map(process_class)).flatten()


def accuracyProcedure(training_data, image_to_classify, predictors, band_names,
                       ancillary, classifier, class_property=None, seed=None, train_prop=None):
    """
    Calculate accuracy metrics using a subset of the training data.
    @param {ee.FeatureCollection} training_data     Training data
    @param {ee.Image}             image_to_classify CCDC coefficient stack to classify
    @param {array}                predictors        List of predictor variables as strings
    @param {array}                band_names        List of band names to classify
    @param {array}                ancillary         List of ancillary predictor data
    @param {ee.Classifier}        classifier        Earth engine classifier with parameters
    @param {string}               class_property    Attribute name with land cover label
    @param {number}               seed              Seed for random column generator
    @param {float}                train_prop        Proportion of data to use for training subset
    @returns {ee.ConfusionMatrix} Confusion matrix as calculated by the train/test subset
    """
    import math as _math
    import random as _random
    train_prop = .4
    seed = seed or _math.ceil(_random.random() * 1000)
    class_property = 'LC_Num'
    training_data = training_data.randomColumn('random', seed).sort('random')
    training_data = subsetTraining(training_data, train_prop, seed, class_property)
    test_subset_test = training_data.filterMetadata('train', 'equals', 0)
    test_subset_train = training_data.filterMetadata('train', 'equals', 1)

    input_list = getInputFeatures(1, image_to_classify, predictors, band_names, ancillary)
    input_features = ee.List(input_list).get(0)

    # Train the classifier
    trained = classifier.train(**{
        'features': test_subset_train,
        'classProperty': class_property,
        'inputProperties': input_features
    })
    classified = test_subset_test.classify(trained)
    conf_matrix = classified.errorMatrix(class_property, 'classification')
    return conf_matrix


def classifyCoefs(image_to_classify, band_names, ancillary, ancillary_features,
                   training_data, classifier, study_area=None, class_property=None,
                   coefs=None, train_prop=None, seed=None):
    """
    Classify single set of CCDC coefficients. Useful for quick parameter testing and debugging.
    @param {ee.Image}             image_to_classify  Single set of ccdc coefficients to classify
    @param {array}                band_names         List of band names to classify
    @param {array}                ancillary          List of ancillary predictor data
    @param {ee.Image}             ancillary_features Ancillary data image
    @param {ee.FeatureCollection} training_data      Training data
    @param {ee.Classifier}        classifier         Earth engine classifier with parameters
    @param {ee.Geometry}          study_area         Boundaries of region to subset training data
    @param {string}               class_property     Attribute name with land cover label
    @param {array}                coefs              List of coefficients to classify
    @param {float}                train_prop         Proportion of data to use for training subset
    @param {number}               seed               Seed for random column generator
    @returns {ee.Image} Classified image
    """
    train_prop = train_prop or None
    study_area = study_area or None
    training_data = ee.FeatureCollection(training_data)
    image_to_classify = ee.Image(image_to_classify)

    # Subset training data to study area if specified
    if study_area:
        training_data = training_data.filterBounds(study_area)

    # Test withholding subset of data and classifying
    if train_prop:
        accuracyProcedure(training_data, seed, train_prop)

    # Predictor names selected for classification.
    predictors = (ee.List(band_names).map(
        lambda b: ee.List(coefs).map(
            lambda i: ee.String(b).cat('_').cat(i)))
        .flatten().cat(ancillary_features))

    # Train the classifier
    trained = classifier.train(**{
        'features': training_data,
        'classProperty': class_property,
        'inputProperties': predictors
    })

    bands = image_to_classify.addBands(ancillary)
    classified = bands.select(predictors).classify(trained).int()
    return classified


def classifySegments(image_to_classify, number_of_segments, band_names, ancillary,
                      ancillary_features, training_data, classifier, study_area=None,
                      class_property=None, coefs=None, train_prop=None, seed=None,
                      subset_training=True):
    """
    Classify stack of CCDC coefficient, band-separated by segment.
    @param {ee.Image}             image_to_classify    CCDC coefficient stack to classify
    @param {number}               number_of_segments   Number of segments to classify
    @param {array}                band_names           List of band names to classify
    @param {array}                ancillary            List of ancillary predictor data
    @param {ee.Image}             ancillary_features   Ancillary data image
    @param {ee.FeatureCollection} training_data        Training data
    @param {ee.Classifier}        classifier           Earth engine classifier with parameters
    @param {ee.Geometry}          study_area           Boundaries of region to subset training data
    @param {string}               class_property       Attribute name with land cover label
    @param {array}                coefs                List of coefficients to classify
    @param {float}                train_prop           Proportion of data to use for training subset
    @param {number}               seed                 Seed for random column generator
    @param {boolean}              subset_training      True to subset training to geometry
    @returns {ee.Image} Classified stack of CCDC segments
    """
    train_prop = train_prop or None
    study_area = study_area or None
    training_data = ee.FeatureCollection(training_data)
    image_to_classify = ee.Image(image_to_classify)

    # Subset training data to study area if specified
    if study_area and subset_training is not False:
        training_data = training_data.filterBounds(study_area)

    # Test withholding subset of data and classifying
    if train_prop:
        accuracyProcedure(training_data, seed, train_prop)

    # Input bands
    predictors = (ee.List(band_names).map(
        lambda b: ee.List(coefs).map(
            lambda i: ee.String(b).cat('_').cat(i)))
        .flatten().cat(ancillary_features))
    input_list = getInputFeatures(1, image_to_classify, predictors, band_names, ancillary)
    input_features = input_list[0]

    # Train the classifier
    trained = classifier.train(**{
        'features': training_data,
        'classProperty': class_property,
        'inputProperties': input_features
    })

    # Map over segments
    def classify_segment(seg):
        seg_input_list = getInputFeatures(seg, image_to_classify, predictors, band_names, ancillary)
        seg_input_features = seg_input_list[0]
        seg_bands = seg_input_list[1]
        seg_str = ee.String('S').cat(ee.String(ee.Number(seg).int8()))
        class_name = seg_str.cat('_classification')
        start_name = seg_str.cat('_tStart')
        return (seg_bands.select(seg_input_features)
                         .classify(trained)
                         .updateMask(image_to_classify.select(start_name).neq(0))
                         .rename([class_name])
                         .int())

    segments_classified = ee.List.sequence(1, number_of_segments).map(classify_segment)

    # segmentsClassified is returned as a list so first convert to Collection
    classified = ee.ImageCollection(segments_classified)

    # When reducing to bands the names change and gives an error upon export
    bns = (ee.List(classified
                   .map(lambda i: i.set('bn', i.bandNames()))
                   .aggregate_array('bn'))
             .flatten())

    # Reduce to bands and rename to original band names
    classified = classified.toBands().rename(bns)
    return classified


def parseConfMatrix(im, attribute='confMatrix'):
    """
    Parse confusion matrix from string.
    @param {ee.Image} im         Classified image with confusion matrix in metadata
    @param {string}   attribute  Name of attribute with confusion matrix
    @returns {ee.Image} Image with confusion matrix parsed and users/producers accuracy set
    """
    # Parse confusion matrix
    conf = ee.String(im.get(attribute))
    conf = conf.slice(1).slice(0, -2)
    split = conf.split('],').map(lambda lst: ee.String(lst).slice(1).split(',').map(
        lambda s: ee.Number.parse(s)))
    conf_matrix = ee.ConfusionMatrix(ee.Array(split))

    # Now create dictionary of users and producers accuracy
    users = conf_matrix.consumersAccuracy().project([1]).toList()
    keys = ee.List.sequence(0, users.length().subtract(1)).map(
        lambda num: ee.String(ee.Number(num).int8()))
    names = keys.map(lambda key: ee.String('users_class_').cat(key))
    users_dict = ee.Dictionary.fromLists(names, users)

    producers = conf_matrix.producersAccuracy().project([0]).toList()
    keys = ee.List.sequence(0, users.length().subtract(1)).map(
        lambda num: ee.String(ee.Number(num).int8()))
    names = keys.map(lambda key: ee.String('producers_class_').cat(key))
    producers_dict = ee.Dictionary.fromLists(names, producers)
    im = im.setMulti(producers_dict.combine(users_dict))
    return im


def getLcAtDate(segs, date, number_of_segments=None, ccdc_image=None, behavior=None,
                band_names=None, input_features=None, spec_image=None, date_format=None):
    """
    Calculate landcover at a date based on pre-classified segments.
    @param {ee.Image} segs                Classified ccd segment image
    @param {string}   date                Date of land cover to retrieve in format 'YYYY-MM-DD'
    @param {number}   number_of_segments  Number of segments in classification image
    @param {string}   ccdc_image          Image with CCDC results
    @param {string}   behavior            Behavior when date is between segments ('none','before','after')
    @param {array}    band_names          List of band names (such as 'BLUE','GREEN')
    @param {array}    input_features      List of input feature names (such as 'INTP' and 'RMSE')
    @param {ee.Image} spec_image          Pre-built CCD image (skips buildCcdImage if provided)
    @param {number}   date_format         Date format of the input date
    @returns {ee.Image} matchingDate  Landcover classification image at date specified in parameter
    """
    segs = ee.Image(segs)
    band_names = band_names or ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP']
    input_features = input_features or ['INTP', 'SLP', 'PHASE', 'AMPLITUDE', 'COS', 'SIN', 'COS2', 'SIN2']
    number_of_segments = number_of_segments or ee.Image(segs).bandNames().length()
    behavior = behavior or 'after'

    if date_format == 0:
        date_format = 0
    elif date_format and date_format > 0:
        date_format = date_format
    else:
        date_format = 1

    # Turn array image into image
    spec_image = spec_image or ccdc_utils.buildCcdImage(ccdc_image, number_of_segments, band_names)
    t_starts = spec_image.select('.*tStart')
    t_ends = spec_image.select('.*tEnd')

    date_formatted = date_utils.convertDate({
        'inputFormat': 3,
        'inputDate': date,
        'outputFormat': date_format
    })

    if behavior == 'before':
        date_mask = t_starts.lt(date_formatted)
        matching_date = segs.updateMask(date_mask).reduce(ee.Reducer.lastNonNull())
    elif behavior == 'after':
        date_mask = t_ends.gt(date_formatted)
        matching_date = segs.updateMask(date_mask).reduce(ee.Reducer.firstNonNull())
    else:
        date_mask = t_starts.lt(date_formatted).And(t_ends.gt(date_formatted))
        matching_date = segs.updateMask(date_mask).reduce(ee.Reducer.firstNonNull())

    return matching_date


def getMode(folder, matching_string):
    """
    Get mode classification from a stack of overlapping result files.
    @param {string} folder           Path to the folder containing the result files
    @param {string} matching_string  Identifier in the result file names
    @returns {ee.Image} Band-wise mode classification
    """
    lst = ee.data.getList({'id': folder})
    ims = []
    for item in lst:
        asset_id = item['id']
        if matching_string in asset_id:
            ims.append(ee.Image(asset_id))
    return ee.ImageCollection(ims).reduce(ee.Reducer.mode())


def getInputDict(band_names, input_features, ancillary_features):
    """
    Build a dictionary of which inputs were used.
    @param {array} band_names         List of spectral band names used
    @param {array} input_features     List of CCDC features used
    @param {array} ancillary_features List of ancillary feature names used
    @returns {list} [inputDict, predictors]
    """
    all_possible_inputs = [
        'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7',
        'BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP', 'INTP',
        'AMPLITUDE', 'PHASE', 'AMPLITUDE_1', 'PHASE_1', 'AMPLITUDE_2',
        'PHASE_2', 'SLP', 'COS', 'SIN', 'COS2', 'SIN2', 'COS3', 'SIN3',
        'RMSE', 'ELEVATION', 'ASPECT', 'DEM_SLOPE', 'RAINFALL', 'TEMPERATURE',
        'AMPLITUDE2', 'PHASE2', 'AMPLITUDE3',
        'PHASE3', 'WATER_OCCURRENCE', 'POPULATION', 'TREE_COVER'
    ]
    all_actual_inputs = band_names + input_features + ancillary_features
    input_booleans = [inp in all_actual_inputs for inp in all_possible_inputs]
    input_dict = ee.Dictionary.fromLists(all_possible_inputs, input_booleans)
    predictors = (ee.List(band_names).map(
        lambda b: ee.List(input_features).map(
            lambda i: ee.String(b).cat('_').cat(i)))
        .flatten().cat(ancillary_features))
    return [input_dict, predictors]


def getPredictors(band_names, input_features, ancillary_features):
    """
    Get list of predictor variable names.
    @param {array} band_names         List of spectral band names
    @param {array} input_features     List of CCDC feature names
    @param {array} ancillary_features List of ancillary feature names
    @returns {ee.List} Flattened list of all predictor names
    """
    return (ee.List(band_names).map(
        lambda b: ee.List(input_features).map(
            lambda i: ee.String(b).cat('_').cat(i)))
        .flatten().cat(ancillary_features))


def loadResults(result_format, change_results, study_region, segs, band_names):
    """
    Load the CCD results as an image with bands corresponding to CCD segments.
    @param {string}       result_format   Format of the results ('SegImage', 'SegCollection', other)
    @param {string}       change_results  Asset path or image collection of results
    @param {ee.Geometry}  study_region    Region to filter collection
    @param {list}         segs            Segment names list
    @param {list}         band_names      Band names list
    @returns {ee.Image} CCD image in long format
    """
    if result_format == 'SegImage':
        ccd_image = ee.Image(change_results)
    elif result_format == 'SegCollection':
        ccd_image = ee.ImageCollection(change_results).filterBounds(study_region).mosaic()
    else:
        coef_image = ee.ImageCollection(change_results).filterBounds(study_region).mosaic()
        ccd_image = ccdc_utils.buildCcdImage(coef_image, len(segs), band_names)
    return ccd_image


def getLC(img, date):
    """
    Get landcover at a specific date, combining before and after lookups.
    @param {ee.Image} img   Classified segment image
    @param {string}   date  Date string in 'YYYY-MM-DD' format
    @returns {ee.Image} Single-band landcover image at the given date
    """
    date_classification_after = getLcAtDate(img, date, None, None, 'after',
                                             None, None, None, 1)
    date_classification_before = getLcAtDate(img, date, None, None, 'before',
                                              None, None, None, 1)
    date_classification = (ee.Image.cat([date_classification_after, date_classification_before])
                              .reduce(ee.Reducer.firstNonNull()))
    return date_classification.rename([ee.String(date)])


def makeYearlyMaps(results, years=None, month=None, day=None):
    """
    Make a list of yearly land cover images.
    @param {ee.Image} results  Classified segment image
    @param {ee.List}  years    List of years to map over
    @param {number}   month    Month for the target date
    @param {number}   day      Day for the target date
    @returns {ee.List} List of land cover images
    """
    years = years or ee.List.sequence(2000, 2020)
    month = month or 6
    day = day or 1

    def format_date(y):
        p1 = ee.String(ee.Number(y).int())
        p2 = ee.String('-')
        p3 = ee.String(ee.Number(month).int())
        p4 = ee.String(ee.Number(day).int())
        return p1.cat(p2).cat(p3).cat(p2).cat(p4)

    formatted = years.map(format_date)

    def make_map(y):
        return getLC(results, y).rename('lc')

    ims = ee.List(formatted).map(make_map)
    return ee.List(ims)


def simpleClassification(fc, atts=None, prop=None, classifier=None):
    """
    Simple classification of CCDC coefficients.
    Allows for classification of the first CCDC segment and only requires
    a training data feature collection.
    @param {ee.FeatureCollection} fc          Training data feature collection
    @param {ee.List}              atts        Attribute names to use
    @param {string}               prop        Class property name
    @param {ee.Classifier}        classifier  Earth engine classifier
    @returns {ee.Image} Classified image
    """
    atts = atts or fc.first().propertyNames().removeAll([
        'ID', 'Start_Year', 'dataPath', 'End_Year', 'Level1_Ecoregion', 'landcover',
        'Dataset', 'system:index', 'LC_Class', 'Continent_Code', 'Glance_Class_ID_level1',
        'Glance_Class_ID_level2', 'Level1_Ecoregion', 'Level2_Ecoregion', 'Dataset_Code',
        'Continent', 'Middle_Year', 'trainYear'])
    classifier = classifier or ee.Classifier.smileRandomForest(200)
    prop = prop or 'Glance_Class_ID_level1'

    ancillary_features = ['ELEVATION', 'ASPECT', 'DEM_SLOPE', 'RAINFALL', 'POPULATION',
                           'WATER_OCCURRENCE']
    ancillary = input_utils.getAncillary()
    band_names = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'TEMP']
    coefs = ['INTP', 'SLP', 'COS', 'SIN', 'RMSE', 'COS2', 'SIN2', 'COS3', 'SIN3']

    ccdc_collection_filtered = (ee.ImageCollection('projects/CCDC/v2')
                                   .filterMetadata('system:index', 'starts_with', 'z_')
                                   .mosaic())
    ccd_image = ccdc_utils.buildCcdImage(ccdc_collection_filtered, 1, band_names)
    predictors = (ee.List(band_names).map(
        lambda b: ee.List(coefs).map(
            lambda i: ee.String(b).cat('_').cat(i)))
        .flatten().cat(ancillary_features))

    input_list = ee.List(getInputFeatures(1, ccd_image, predictors, band_names, ancillary))
    input_features = ee.List(input_list.get(0))
    image_to_classify = ee.Image(input_list.get(1))

    # Train the classifier
    trained = classifier.train(**{
        'features': fc,
        'classProperty': prop,
        'inputProperties': input_features
    })
    return image_to_classify.classify(trained)
