import sys
import os
sys.path.append('%s/data/python' % os.environ['ANTELOPE'])
from antelope.datascope import closing, dbopen
from antelope.stock import str2epoch, epoch2str
from misc_tools import create_event_list
from eqloc3d_classes import *
#with closing(dbopen('/Users/mcwhite/staging/dbs/anza_sub/anza', 'r')) as db:
#    tbl_event = db.schema_tables['event']
#    view1 = tbl_event.join('origin')
#    view1 = view1.subset('time >= _2013315 00:00:00_')
#    print view1.record_count
#    view1 = view1.separate('event')
#    event_list = create_event_list(view1, 'CSS3.0')
#    origin = event_list[5].preferred_origin
#    print origin.lat

event_list = create_event_list('pha.phase_210to9_july2011.dat', 'SCEDC')
len(event_list)
print event_list[0]
