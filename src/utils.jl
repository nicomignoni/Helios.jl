# Bounds for latitude and longitude
const LATITUDE_MAX = 90.0
const LATITUDE_MIN = -90.0
const LONGITUDE_MAX = 180.0
const LONGITUDE_MIN = -180.0

const LATITUDE_RANGE = LATITUDE_MAX - LATITUDE_MIN
const LONGITUDE_RANGE = LONGITUDE_MAX - LONGITUDE_MIN

function interpolated_value(matrix, latitude, longitude, I...)
    nrows, ncols = size(matrix)
    lat_ratio = nrows * (latitude - LATITUDE_MIN) / LATITUDE_RANGE
    lon_ratio = ncols * (longitude - LONGITUDE_MIN) / LONGITUDE_RANGE

    lat_index_up, lat_index_down = Int(ceil(lat_ratio)), Int(floor(lat_ratio))
    lon_index_up, lon_index_down = Int(ceil(lon_ratio)), Int(floor(lon_ratio))
    
    # Index position between rounded row and column as [0, 1] fraction
    lat_frac, lon_frac = lat_ratio - lat_index_down, lon_ratio - lon_index_down

    # Bilinear interpolation of row-column data
    return lat_frac * lon_frac * matrix[lat_index_down, lon_index_down, I...] + 
           lat_frac * (1 - lon_frac) * matrix[lat_index_down, lon_index_up, I...] + 
           (1 - lat_frac) * lon_frac * matrix[lat_index_up, lon_index_down, I...] + 
           (1 - lat_frac) * (1 - lon_frac) * matrix[lat_index_up, lon_index_up, I...] 
end
