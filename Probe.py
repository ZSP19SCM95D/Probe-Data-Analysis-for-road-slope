from datetime import datetime, date, time
from collections import defaultdict
import math
from scipy.stats import norm
import sys


class ReadFile:
    def __init__(self, filedir):
        self.dir = filedir
        self.file = open(filedir, 'r', encoding='utf8')

    def readLinks(self):

        mp = GerMap()

        line = self.file.readline()
        while len(line):

            ftrs = line.split(',')
            latlonst1 = [ll.split('/') for ll in ftrs[14].split('|')]
            latlons = list()
            for ll in latlonst1:
                sl0 = float(ll[0])
                sl1 = float(ll[1])
                if len(ll[2]):
                    sl2 = float(ll[2])
                else:
                    sl2 = 0
                latlons.append((sl0, sl1, sl2))
            rnode = Node(int(ftrs[1]), latlons[0][0], latlons[0][1], latlons[0][2])
            nrnode = Node(int(ftrs[2]), latlons[-1][0], latlons[-1][1], latlons[-1][2])
            link = Link(int(ftrs[0]), rnode, nrnode, float(ftrs[3]))
            for i in range(1, len(latlons) - 1):
                link.setShape(latlons[i][0], latlons[i][1], latlons[i][2])

            if len(ftrs[15]):
                dstcurt1 = [ll.split('/') for ll in ftrs[15].split('|')]
                dstcurs = list()
                for ll in dstcurt1:
                    sl0 = float(ll[0])
                    sl1 = float(ll[1])
                    dstcurs.append((sl0, sl1))
                for d, c in dstcurs:
                    link.setCurvature(d, c)

            if len(ftrs[16]) > 1:
                dstslpt1 = [ll.split('/') for ll in ftrs[16].split('|')]
                dstslps = list()
                for ll in dstslpt1:
                    sl0 = float(ll[0])
                    sl1 = float(ll[1])
                    dstslps.append((sl0, sl1))
                for d, s in dstslps:
                    link.setSlope(d, s)
            mp.addLink(link)

            line = self.file.readline()

        return mp

    def readProbes(self, maxlength=5000):
        prbs = defaultdict(list)
        for i in range(maxlength):
            line = self.file.readline()
            if line:
                ftrs = line.split(',')
                sid = int(ftrs[0])
                dtt1 = ftrs[1].split()
                dt2 = dtt1[0].split('/')
                d = date(int(dt2[2]), int(dt2[0]), int(dt2[1]))

                tt3 = dtt1[1].split(':')
                if dtt1[-1] == 'AM' or int(tt3[0]) == 12:
                    t = time(int(tt3[0]), int(tt3[1]), int(tt3[2]))
                elif dtt1[-1] == 'PM':
                    t = time(int(tt3[0]) + 12, int(tt3[1]), int(tt3[2]))
                else:
                    continue
                dt = datetime.combine(d, t)

                prb = Probe(sid, dt, int(ftrs[2]), float(ftrs[3]), float(ftrs[4]), float(ftrs[5]), float(ftrs[6]),
                            float(ftrs[7]))
                prbs[sid].append(prb)
            else:
                return prbs

        return prbs

class GerMap:
	def __init__(self):
		self.ERADIUS = 6371000
		self.MAXLAT = 53.4375
		self.MINLAT = 50.6250
		self.MAXLON = 11.25
		self.MINLON = 8.4375
		self.UNIT = 0.003
		self.index = defaultdict(list)
		self.lincolle = dict()

	def calXY(self, lat, lon):
		x = int((lon - self.MINLON) / self.UNIT)
		y = int((lat - self.MINLAT) / self.UNIT)

		return x, y


	def addLink(self, link):
		self.lincolle[link.id] = link
		self.addLinkIdx(link)

	def addLinkIdx(self, link):
		lat, lon = link.MidPoint()
		x , y = self.calXY(lat, lon)
		self.index[(x, y)].append(link.id)

	def searchLinks(self, lat, lon):
		x, y = self.calXY(lat, lon)
		lids = self.index[(x, y)]
		return [self.lincolle[id] for id in lids]

class Link:
	def __init__(self, lid, rnode, nrnode, leng):
		self.id = lid
		self.refNode = rnode
		self.nrefNode = nrnode
		self.length = leng
		self.shapeInfo = list()
		self.shapeInfo.append((rnode.lat, rnode.lon, rnode.elev))
		self.shapeInfo.append((nrnode.lat, nrnode.lon, nrnode.elev))
		self.curvatureInfo = None
		self.slopeInfo = None

	def setShape(self, lat, lon, elev = None):
		self.shapeInfo.insert(-1, (lat, lon, elev))

	def setCurvature(self, dist, curvature):
		if not self.curvatureInfo:
			self.curvatureInfo = list()
		self.curvatureInfo.append((dist, curvature))

	def setSlope(self, dist, slope):
		if not self.slopeInfo:
			self.slopeInfo = list()
		self.slopeInfo.append((dist, slope))

	def getElev(self, lat, lon):
		for i in range(len(self.shapeInfo)):
			si = self.shapeInfo[i]
			if lat == si[0] and lon == si[1]:
				return si[2], i
		return 0, -1



	def calSegProjection(self, lat, lon, lat0, lon0, lat1, lon1):
		y0 = lat0
		x0 = lon0
		y1 = lat1 - y0
		x1 = math.cos(lat1/ 180 * math.pi) * (lon1 - x0)
		y2 = lat - y0
		x2 = math.cos(lat / 180 * math.pi) * (lon - x0)
		x = (x1 * x2 + y1 * y2) / (x1 + y1)
		if x1:
			y = y1 * x / x1
		else:
			y = y1 * x / 0.0001
		lonx = x / math.cos(y / 180 * math.pi)
		if x < 0:
			d = Node.dist(lat, lon, lat0, lon0)
		elif x > x1:
			d = Node.dist(lat, lon, lat1, lon1)
		else:
			d = Node.dist(y, lonx, lat, lon)
		return y + lat0, lonx + lon0, d, self.id, lat0, lon0

	def calProjDist(self, lat, lon):
		tem = list()
		for i in range(len(self.shapeInfo) - 1):
			p = self.calSegProjection(lat, lon, self.shapeInfo[i][0], self.shapeInfo[i][1], self.shapeInfo[i + 1][0], self.shapeInfo[i + 1][1])
			tem.append(p)
		return sorted(tem, key = lambda s:s[2])

	def MidPoint(self):
		return (self.shapeInfo[0][0] + self.shapeInfo[-1][0]) / 2 , (self.shapeInfo[0][1] + self.shapeInfo[-1][1]) / 2

class Match:
	def __init__(self, plist, mp):
		self.prbls = plist
		self.map = mp

	def probeCandidates(self):
		tre = list()
		for pb in self.prbls:
			lks = self.map.searchLinks(pb.lat, pb.lon)
			cps = list()
			for lk in lks:
				cpst1 = lk.calProjDist(pb.lat, pb.lon)
				cps.extend(cpst1)
			tre.append(cps)

		return tre

	def Vc0c1(self, p0, p1, cp0, cp1):
		nmn = Node.dist(p0.lat, p0.lon, p1.lat, p1.lon)
		dnm = Node.dist(cp0[0], cp0[1], cp1[0], cp1[1])
		if dnm:
			return nmn / dnm
		else:
			return nmn / 0.01



	def findMatchedSequence(self, pcps):
		f = defaultdict(float)
		pre = dict()
		for cp in pcps[0]:
			f[cp] = norm.pdf(cp[2], 0, 20)
		for i in range(1, len(pcps)):
			for ccp in pcps[i]:
				print('\t%f,%f,%f,%d,%f,%f' % ccp)
				mx = 0
				for prp in pcps[i - 1]:
					alt = f[prp] + norm.pdf(ccp[2], 0, 20) * self.Vc0c1(self.prbls[i - 1], self.prbls[i], list(prp), list(ccp))
					if alt > mx:
						mx = alt
						pre[ccp] = prp
					f[ccp] = mx

		ptat1 = sorted([(lp, f[lp]) for lp in pcps[-1]], key = lambda x:-x[1])
		if not len(ptat1):
			return list()
		pta = ptat1[0][0]
		re = list()
		for i in range(len(pcps) - 1, 0, -1):

			re.insert(0, pta)

			if pta in pre.keys():
				pta = pre[pta]
			else:
				return zip(self.prbls[-len(re):], re)

		re.insert(0, pta)

		return zip(self.prbls, re)

class Node:
	def __init__(self, nid, latitude, longitude, elevation = None):
		self.id = nid
		self.lat = latitude
		self.lon = longitude
		self.elev = elevation
		pass

	@staticmethod
	def dist(lat1, lon1, lat2, lon2):
		return math.asin((math.sin((lat2 - lat1) / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin((lat2 - lat1) / 2) ** 2) ** 0.5)

	@staticmethod
	def degree(lat0, lon0, lat1, lon1):
		y = lat1 - lat0
		if not y:
			return 90.0
		x = math.cos(lat1 / 180 * math.pi) * (lon1 - lon0)
		return math.atan(x / y) / math.pi * 180

	@staticmethod
	def slope(lat0, lon0, elev0, lat1, lon1, elev1):
		return 180 / math.pi * math.atan(Node.dist(lat0, lon0, lat1, lon1) / (elev1 - elev0))

class Probe:
	def __init__(self, sid , dtime, scode, latitude, longitude, altitude, spd, hding):
		self.id = sid
		self.sourceCode = scode
		self.dttm = dtime
		self.lat = latitude
		self.lon = longitude
		self.alt = altitude
		self.speed = spd
		self.heading = hding

class Slope:
	def __init__(self, pbcd, mp):
		self.probcan = pbcd
		self.map = mp

	def calSlope(self):
		re = list()
		for i in range(len(self.probcan)):
			pd = self.probcan[i]
			lid = pd[1][3]
			lk = self.map.lincolle[lid]
			rlat = pd[1][4]
			rlon = pd[1][5]
			if lk.slopeInfo:

				relev , j = lk.getElev(rlat, rlon)
				slc = Node.slope(rlat, rlon, relev, pd[0].lat, pd[0].lon, pd[0].alt)
				sl = lk.slopeInfo[j][1]
				re.append((rlat, rlon, sl, slc))
			else:
				re.append((rlat, rlon, 0.0, 0.0))

		return re

	def evalSlope(self, slslc):
		return [sl - slc for sl, slc in slslc]

class WriteFile:
	def __init__(self):
		self.DIRMP = 'Partition6467MatchedPoints.csv'
		self.DIRS = 'Partition6467Slopes.csv'

	def writeMP(self, mtsq):
		file = open(self.DIRMP, 'w', encoding = 'utf8')
		for mq in mtsq:
			hd = mq[0].heading
			dg = Node.degree(mq[1][4], mq[1][5], mq[1][0], mq[1][1])
			if abs(hd - dg) > 90:
				direction = 'T'
			else:
				direction = 'F'
			dfr = Node.dist(mq[1][4], mq[1][5], mq[1][0], mq[1][1])
			file.write('%d,%s,%d,%f,%f,%f,%f,%f,%d,%s,%f,%f\n' % (mq[0].id, mq[0].dttm, mq[0].sourceCode, mq[0].lat, mq[0].lon, mq[0].alt, mq[0].speed, mq[0].heading, mq[1][3], direction, mq[1][2], dfr))
		file.close()


	def writeS(self, slp):
		file = open(self.DIRS, 'w', encoding = 'utf8')
		for s in slp:
			file.write('%f,%f,%f,%f\n' % s)
		file.close()


def main():
    readlink = ReadFile('Partition6467LinkData.csv')
    mp = readlink.readLinks()
    readprobe = ReadFile('Partition6467ProbePoints.csv')

    if (len(sys.argv) > 1):
        print('Start reading %d probes ...' % int(sys.argv[1]))

        prbs = readprobe.readProbes(int(sys.argv[1]))
    else:
        print('Start reading %d probes ...' % 5000)

        prbs = readprobe.readProbes()

    sampleIds = prbs.keys()
    mtsq = list()
    for sid in sampleIds:
        mc = Match(prbs[sid], mp)
        pcps = mc.probeCandidates()
        mtsq.extend(mc.findMatchedSequence(pcps))

    wf = WriteFile()

    wf.writeMP(mtsq)
    slope = Slope(mtsq, mp)
    slp = slope.calSlope()
    wf.writeS(slp)
    print('Slopes writing complete ...')


if __name__ == '__main__':
    main()
