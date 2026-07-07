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

Miscellaneous utility functions for working with GEE assets.
"""

import ee


# =========================
# 1. FUNCTIONS ==============
# =========================

def assetsToCollection(folder_path, asset_type, string_match):
    """
    Return a list of assets from a folder that match a given string.
    Filters in ee.data.listAssets are only applied to ImageCollection assets
    and are ignored for Folder assets.

    @param {String} folder_path   Path to the asset folder
    @param {String} asset_type    Type of asset: 'featureCollection', 'Image', or 'object'
    @param {String} string_match  Substring to match against asset IDs
    @returns {list} List of ee objects or IDs matching the filter
    """
    asset_list = ee.data.listAssets(folder_path, {'view': 'BASIC'}).get('assets', [])
    ims = []
    if not asset_list:
        asset_list = []
    for item in asset_list:
        asset_id = item['id']
        if string_match in asset_id:
            if asset_type == 'featureCollection':
                im = ee.FeatureCollection(asset_id).set('assetId', asset_id)
            elif asset_type == 'Image':
                im = ee.Image(asset_id)
            elif asset_type == 'object':
                im = asset_id
            else:
                im = asset_id
            ims.append(im)
    return ims
