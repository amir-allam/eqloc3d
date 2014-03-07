from misc_tools import _Event

event = _Event('~/staging/dbs/anza_sub/anza', 83400)
print event
print  '!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print event.lat, event.lon
print  '!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print event.arrivals
print  '!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print event.arrivals[2].sta, event.arrivals[2].phase
