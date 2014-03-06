from antelope.datascope import closing, dbopen
from misc_tools import StationList, process_events

def _main():
    args = _parse_args()
    params = _parse_pf(args.pf)
    while True:
        last_lddate = 0
        events = []
        station_list = StationList(dbin)
        with closing(dbopen(args.dbin, 'r')) as dbin:
            tbl_event = dbin.schema_table['event']
            tbl_event = tbl_event.subset('lddate > _%f_' % last_lddate)
            for record in tbl_event.iter_record():
                events.append(record.getv('evid')[0])
        process_events(dbin, dbout, events, station_list)
    return 0

def _parse_args():
    """
    Parse command line options.
    Return dictionary-like object containing arguments.
    """
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('dbin', type=str, help='Input databse.')
    parser.add_argument('dbout', type=str, help='Output databse.')
    parser.add_argument('-p', '--pf', type=str, help='Parameter file.')
    return parser.parse_args()

def _parse_pf(pf):
    """
    Parse parameter file.
    Return dictionary-like object containing parameters.
    """
    from antelope.stock import pfread, pfin
    pf = pfin(pf) if pf else pfread(sys.argv[0])
    return _recursive_eval(pf.pf2dict())

def _recursive_eval(a_dict):
    for key in a_dict:
        if isinstance(a_dict[key], dict):
            a_dict[key] = _recursive_eval(a_dict[key])
        else:
            try:
                a_dict[key] = eval(a_dict[key])
            except Exception:
                pass
        return a_dict

if __name__ == '__main__': sys.exit(_main())
else:
    print '%s - Not a module to import!!' % sys.argv[0]
    sys.exit(-1)
