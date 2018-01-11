import math
import numpy as np

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

class Source(object):
	def __init__(self, name="ERROR", RA=0.0, DEC=0.0):
		self.name = name
		self.RA = RA
		self.DEC = DEC
	def __str__(self):
		return \
			"Source " + self.name + " " +\
			"RA: " + str(self.RA) + " " +\
			"DEC: " + str(self.DEC)
class Event(object):
	def __init__(self, ID=-1, RA=0.0, DEC=0.0, angular_error=0.0):
		self.ID = ID
		self.RA = RA
		self.DEC = DEC
		self.angular_error = angular_error
	def __str__(self):
		return "Event " + str(self.ID) + " " + \
		"RA: " + str(self.RA) + " " + \
		"DEC: " + str(self.DEC) + " " + \
		"Angular Error: " + str(self.angular_error)

sources = {}
#read in the sources
filename = "data/HAWCcatalog.txt" #str(raw_input("Enter Filename >"))
with open(filename, "r") as file:
	for line in file.readlines():
		raw = line.strip("\n").replace("*", "").replace("\xe2\x97\xa6", "").split(" ")
		name = raw[0] + " " + raw[1]
		#convert to radians on readin
		ra = float(raw[4]) * (math.pi / 180.0)
		dec = float(raw[5]) * (math.pi / 180.0)
		sources[name] = Source(name, ra, dec)

"""print the sources
for name in sources:
	print sources[name]
#"""
events = {}
#read in the events
filename = "data/eventsummary_4years.txt"
with open(filename, "r") as file:
	for line in file.readlines():
		if "#" in line: continue
		raw = line.split()
		eventID = int(raw[0])
		dec = float(raw[5]) * (math.pi / 180.0)
		ra = float(raw[6]) * (math.pi / 180.0)
		ang_error = float(raw[7]) * (math.pi / 180.0)
		events[eventID] = Event(eventID, ra, dec, ang_error)
"""print the events
for eventID in events:
	print events[eventID]
#"""

# step 3 - calculate the angular probablitiy between one source and one neutrino event
# using the event_angular_distribution function below

def sph_dot(th1,th2,phi1,phi2):
	return np.sin(th1)*np.sin(th2)*np.cos(phi1-phi2) + np.cos(th1)*np.cos(th2)

# Implementation of the kent distribution
def event_angular_distribution(event,src):
	#print (event)
	#print (src)
	kappa = 1.0 / (event.angular_error)**2
	log_dist = np.log(kappa) - np.log(2*np.pi) - kappa + kappa*sph_dot(np.pi/2-event.DEC, np.pi/2-src.DEC, event.RA, src.RA)
	return np.exp(log_dist)
#print event_angular_distribution(events[1], sources["2HWC J0631+169"])

# step 4 - calculate angular probability between one source and all neutrino events as a list or array
angular_probability = []
for eventID in events:
	prob = event_angular_distribution(events[eventID],sources["2HWC J0631+169"])
	angular_probability.append(prob)
#print(angular_probability)
def H(source,ns):
	sum = 0.
	N   = len(events)
	for key in events:
		#print(events[key], source)
		sum += np.log((float(ns)/N)*event_angular_distribution(events[key], source)+(1 - float(ns)/N)/(4*np.pi))
	return 2*sum

def calcTS(source,ns):
	return H(source, ns)-H(source,0)

def Hcat(ns):
	sum = 0.0
	M = len(sources)
	N = len(events)
	total = 0.0
	for ID in events:
		total = 0.0
		for src_name in sources:
			total += event_angular_distribution(events[ID], sources[src_name])
		sum += np.log( (float(ns)/N)*total/M + (1.0 - float(ns)/N) / (4.0*np.pi) )
	return 2.0 * sum
#calculate Hcat(0) once
Hcat0 = Hcat(0)
def TScat(ns):
	#replaced the Hcat(0) with a constant value
	return Hcat(ns) - Hcat0
ts = []
print "Starting"
for ns in np.arange(0, 50, 0.1):
	ts.append(TScat(ns))
	print ns, ts[-1]
print "max=", max(ts)