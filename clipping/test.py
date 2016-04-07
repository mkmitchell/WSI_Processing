from ocgis import OcgOperations, RequestDataset
from ocgis.test.base import TestBase


class Test(TestBase):

def test_time_region(self):
    uri = 'C:/testclip/WSI_OCGIS_abdu.1979.nc'
    shp = 'C:/testclip/state.shp'
    rd = RequestDataset(uri=uri)
    calc = [{'func': 'sum', 'name': 'sum'}]
    ops_one = OcgOperations(dataset=rd, output_format='numpy', time_region={'month': [1]},
                                 spatial_operation='clip', geom=shp, calc=calc, calc_raw=True, aggregate=True,
                                 calc_grouping='day', prefix='calc', geom_select_sql_where='STATE_NAME="Alabama"')
    ret_one_month = ops_one.execute()
    ops_two = OcgOperations(dataset=rd, output_format='numpy', time_region={'month': [2]},
                                 spatial_operation='clip', geom=shp, calc=calc, calc_raw=True, aggregate=True,
                                 calc_grouping='day', prefix='calc', geom_select_sql_where='STATE_NAME="Alabama"')
    ret_two_month = ops_two.execute()
    ops_original = OcgOperations(dataset=rd, output_format='numpy', time_region={'month': [1, 2]},
                                 spatial_operation='clip', geom=shp, calc=calc, calc_raw=True, aggregate=True,
                                 calc_grouping='day', prefix='calc', geom_select_sql_where='STATE_NAME="Alabama"')
    ret_original = ops_original.execute()
    desired = ret_original[1]['forcalc'].variables['sum'].value  # 11.580645161290322
    ops_no_time_region = OcgOperations(dataset=rd, output_format='numpy',
                                 spatial_operation='clip', geom=shp, calc=calc, calc_raw=True, aggregate=True,
                                 calc_grouping='day', prefix='calc', geom_select_sql_where='STATE_NAME="Alabama"')
    ret_no_time_region = ops_no_time_region.execute()
    field = ret_no_time_region[1]['forcalc']
    indices = []
    for idx in range(field.temporal.shape[0]):
        the_time = field.temporal.value_datetime[idx]
        if the_time.month in [1, 2]:
            indices.append(idx)
    var_sub = field.variables['sum'][:, indices, :, :, :]
    actual = var_sub.value
    self.assertNumpyAll(actual, desired)
