import pandas as pd
import math
import copy

# Global variables
EarthRad_km = 6371.0

class County:
    def __init__(self,name,boundaryPoints):
        self.name = name
        
        # Need at least three boundary points to make a County
        # Each point should contain a lat and long (2 floats)
        # The first and last points are the same, so have to subtract 1
        # from the count of points provided to get total number of coords
        # that make up this county's border
        CoordError = False
        for coord in boundaryPoints:
            if len(coord) != 2:
                CoordError = True
                break
            elif type(coord[0]) is not float or type(coord[1]) is not float:
                CoordError = True
                break
        if CoordError:
            self.borders = None
            print("Error: There is an issue with the data of the coordinates supplied.")
            return
        elif len(boundaryPoints) < 4:
            self.borders = None
            print("Error: There are not sufficient points to make a county. Need at least three distinct points.")
            return
        elif boundaryPoints[0] != boundaryPoints[-1]:
            self.borders = None
            # checking what point vals are
            print("Boundary point 0: "+str(boundaryPoints[0]))
            print("Boundary point -1: "+str(boundaryPoints[-1]))
            print("Error: Our method for generating a County object requires that the first and last coordinates are the same.")
            return
        # Our data appears to be good, so proceed with creating self.borders
        allBorders = []
        # initialize ranges to be the lat and long values of the first point
        rangeLats = [boundaryPoints[0][0],boundaryPoints[0][0]]
        rangeLongs = [boundaryPoints[0][1],boundaryPoints[0][1]]
        for idx in range(len(boundaryPoints)-1):
            # the idx of the connecting coordinate to this coord is always idx+1
            # note that the last coord is the same as the first
            nextIdx = idx+1
            # check if this point's lat or longs are extremes
            if boundaryPoints[idx][0] > rangeLats[1]:
                rangeLats[1] = boundaryPoints[idx][0]
            elif boundaryPoints[idx][0] < rangeLats[0]:
                rangeLats[0] = boundaryPoints[idx][0]
            if boundaryPoints[idx][1] > rangeLongs[1]:
                rangeLongs[1] = boundaryPoints[idx][1]
            elif boundaryPoints[idx][1] < rangeLongs[0]:
                rangeLongs[0] = boundaryPoints[idx][1]

            # now we create a Segment from these two coordinates
            aBorder = Segment(boundaryPoints[idx],boundaryPoints[nextIdx])
            # and append it to our running list of Segments
            allBorders.append(aBorder)
        # set self.borders to the list of all the borders created by the coord Segments
        self.borders = allBorders
        # set this county's Lat and Long ranges
        self.maxLat = rangeLats[1]
        self.minLat = rangeLats[0]
        self.maxLong = rangeLongs[1]
        self.minLong = rangeLongs[0]

        # Create variables to store COM and Area once computed
        self.COM = [None,None]
        self.Area = None

    # Class method to print all of the borders for this County
    def PrintMyBorders(self):
        s = ''
        s += self.name + ' county has ' + str(len(self.borders)) + ' border segments\n'
        s += 80*'-'
        ctr = 1
        for border in self.borders:
            s += '\nBorder ' + str(ctr) + ':\t'
            s += border.PrintSegmentPoints()
            ctr += 1
        s += '\n'+(80*'-')
        return s

    # Class method to print COM and Area
    def StringCOMandArea(self):
        s = '\n'
        #s += (80*'-')
        #s += '\n'
        s += self.name + " county's Center of Mass is located at " + str(self.COM)
        s += "\n"
        s += self.name + " county has an approximate surface area of " + str(self.Area) + " km^2"
        s += "\n"
        s += (80*'-')
        s += "\n"
        return s

    def GetBorderIdxsSpanningLat(self,aLat):
        outLs = []
        for idx in range(len(self.borders)):
            if self.borders[idx].CrossesLat(aLat):
                outLs.append(idx)
        return outLs
    
    def GetBorderIdxsSpanningLong(self,aLong):
        outLs = []
        for idx in range(len(self.borders)):
            if self.borders[idx].CrossesLong(aLong):
                outLs.append(idx)
        return outLs

    def VertexOnlyTouchesLat(self,aLat,border1,border2,stepAmt):
        if border1.CrossesLat(aLat+stepAmt) and border2.CrossesLat(aLat+stepAmt):
            return True
        elif border1.CrossesLat(aLat-stepAmt) and border2.CrossesLat(aLat-stepAmt):
            return True
        else:
            return False

    def VertexOnlyTouchesLong(self,aLong,border1,border2,stepAmt):
        if border1.CrossesLong(aLong+stepAmt) and border2.CrossesLong(aLong+stepAmt):
            return True
        elif border1.CrossesLong(aLong-stepAmt) and border2.CrossesLong(aLong-stepAmt):
            return True
        else:
            return False

    def ComputeCOMLat(self, stepAmt):
        # create return variable
        retObj = {}
        retObj['COM_Lat'] = None
        retObj['CountyArea'] = None

        lsLats = []
        lsDists = []
        startLat = GetStartingValue(self.maxLat,stepAmt)
        # Progress measuring:
        print("\n\nComputing COM Lat...")
        numIntervals = 50
        numLatChecks = int((self.maxLat - self.minLat)/stepAmt)
        interval = math.floor(numLatChecks/numIntervals)
        intervals = [int(interval*intervalIdx) for intervalIdx in range(numIntervals+1)]
        latProgressCtr = 0
        # Print what full progress bar looks like first for comparison
        print("0% " + ((numIntervals+1)*"*") + " 100% (for reference)")
        print("0% ",end="")
        # while loop where all the magic happens
        curLat = startLat
        while curLat >= self.minLat:
            # Progress printing:
            if latProgressCtr in intervals:
                print("*",end="")
            latProgressCtr += 1
            # identify/extract the border segments that span this lat
            bIdxsAtLat = self.GetBorderIdxsSpanningLat(curLat)
            numBordersFound = len(bIdxsAtLat)
            # get an ordered list of all the Long values for borders that span this Lat
            lsOrderedLongs = []
            for idxOfIdxs in range(numBordersFound):
                # get the longitude value for this border at the curLat
                thisBordIdx = bIdxsAtLat[idxOfIdxs]
                aLong = self.borders[thisBordIdx].getLongAtLat(curLat)
                # check for special case between this border and the prev
                # (as long as this isn't the first border in list and the ordered list of
                # Longs has a value already stored to compare to)
                if idxOfIdxs > 0 and len(lsOrderedLongs) > 0:
                    if lsOrderedLongs[-1] == aLong:
                        prevBordIdx = bIdxsAtLat[idxOfIdxs] - 1
                        if self.VertexOnlyTouchesLat(curLat,self.borders[prevBordIdx],self.borders[thisBordIdx],stepAmt):
                            # if these two borders only just meet up at this latitude, then
                            # we can ignore both of them for purposes of computing the distance
                            # between where they intersect
                            lsOrderedLongs.pop(-1)
                            continue
                        else:
                            # if these two borders do not touch at the given lat, then they must span across it
                            # in which case we accept the one already stored and continue without appending this one
                            continue
                lsOrderedLongs.append(aLong)
            lsOrderedLongs.sort()
            distInside = 0
            if len(lsOrderedLongs) % 2 == 0:
                for l_idx in range(0,len(lsOrderedLongs),2):
                    distInside += (lsOrderedLongs[l_idx+1] - lsOrderedLongs[l_idx])
            if len(lsOrderedLongs) % 2 > 0:
                print("\n\nSomething went wrong because we got "+str(len(lsOrderedLongs))+" crossing longs at "+str(curLat)+" deg Latitude")
            # at end of while statement, decrement curLat by stepAmt
            # and store distInside and Lat
            lsLats.append(curLat)
            lsDists.append(distInside)
            curLat -= stepAmt

        # Progress bar completed:
        print(" 100%")
        ## Prints for visualizing results
        #print("\nPrinting the first 1000 distances inside\n")
        #if len(lsDists) < 1000:
        #    printRng = len(lsDists)
        #else:
        #    printRng = 1000
        #for idx in range(printRng):
        #    print(str(lsDists[idx]))

        # **We need to obtain a third list containing the weights assigned at each lat.
        # The weight assigned is cos(lat)*dist so that at the equator the weight is maximal (1)*dist
        # and at either pole it is minimal (0)*dist
        lsWeights = [math.cos(math.radians(lsLats[i])) * lsDists[i] for i in range(len(lsLats))]
        totWeight = 0
        sumLatWeightProduct = 0
        # perform the summation of the product of lats with their associated weights
        # as well as the total weight to use in the divisor
        for i in range(len(lsWeights)):
            sumLatWeightProduct += lsLats[i] * lsWeights[i]
            totWeight += lsWeights[i]
        # last step is to divide the sum of the Lat-Weight product by the tot Weight.
        # Check to make sure we don't div by 0
        if totWeight != 0:
            retObj['COM_Lat'] = sumLatWeightProduct / totWeight
        # Actually, we need to calculate the area of the county in km^2 so that when
        # calculating the State's COM, we can multiply this county's COM by its appropriate weight.
        # To get this, we multiply our totWeight by the product of the earth's radius squared (in km^2)
        # and the stepAmount (in degrees).
        # Note: we needed to convert from deg to radians in order for an accurate measure. Since
        # our calculation uses two measures in degrees, we need to square this conversion factor.
        retObj['CountyArea'] = totWeight * EarthRad_km**2 * stepAmt * (math.pi/180)**2
        print("\n\nEstimated area for "+self.name+" county is "+str(retObj['CountyArea'])+" km^2")
        return retObj
    
    def ComputeCOMLong(self, stepAmt):
        # create return variable
        retObj = {}
        retObj['COM_Long'] = None
        retObj['CountyArea'] = None

        lsLongs = []
        lsWeights = []
        startLong = GetStartingValue(self.maxLong,stepAmt)
        # Progress measuring:
        print("\n\nComputing COM Long...")
        numIntervals = 50
        numLongChecks = int((self.maxLong - self.minLong)/stepAmt)
        interval = math.floor(numLongChecks/numIntervals)
        intervals = [int(interval*intervalIdx) for intervalIdx in range(numIntervals+1)]
        longProgressCtr = 0
        # Print what full progress bar looks like first for comparison
        print("0% " + ((numIntervals+1)*"*") + " 100% (for reference)")
        print("0% ",end="")
        # while loop where all the magic happens
        curLong = startLong
        while curLong >= self.minLong:
            # Progress printing:
            if longProgressCtr in intervals:
                print("*",end="")
            longProgressCtr += 1
            # identify/extract the border segments that span this lat
            bIdxsAtLong = self.GetBorderIdxsSpanningLong(curLong)
            numBordersFound = len(bIdxsAtLong)
            # get an ordered list of all the Lat values for borders that span this Lat
            lsOrderedLats = []
            for idxOfIdxs in range(numBordersFound):
                # get the latitude value for this border at the curLong
                thisBordIdx = bIdxsAtLong[idxOfIdxs]
                aLat = self.borders[thisBordIdx].getLatAtLong(curLong)
                # check for special case between this border and the prev
                # (as long as this isn't the first border in list and the ordered list of
                # Lats has a value already stored to compare to)
                if idxOfIdxs > 0 and len(lsOrderedLats) > 0:
                    if lsOrderedLats[-1] == aLat:
                        prevBordIdx = bIdxsAtLong[idxOfIdxs] - 1
                        if self.VertexOnlyTouchesLong(curLong,self.borders[prevBordIdx],self.borders[thisBordIdx],stepAmt):
                            # if these two borders only just meet up at this latitude, then
                            # we can ignore both of them for purposes of computing the distance
                            # between where they intersect
                            lsOrderedLats.pop(-1)
                            continue
                        else:
                            # if these two borders do not touch at the given lat, then they must span across it
                            # in which case we accept the one already stored and continue without appending this one
                            continue
                lsOrderedLats.append(aLat)
            lsOrderedLats.sort()
            distInside = 0
            if len(lsOrderedLats) % 2 == 0:
                for l_idx in range(0,len(lsOrderedLats),2):
                    # formula to calculate SA of sphere using integration requires
                    # us to take the difference of the sine of the two lat values
                    distInside += (math.sin(math.radians(lsOrderedLats[l_idx+1])) - math.sin(math.radians(lsOrderedLats[l_idx])))
            if len(lsOrderedLats) % 2 > 0:
                print("\n\nSomething went wrong because we got "+str(len(lsOrderedLats))+" crossing lats at "+str(curLong)+" deg Longitude")
            # at end of while statement, decrement curLong by stepAmt
            # and store distInside and Lat
            lsLongs.append(curLong)
            lsWeights.append(distInside)
            
            # at end of while statement, decrement curLong by stepAmt
            curLong -= stepAmt
        
        # Progress bar completed:
        print(" 100%")
        ## Prints for visualizing results
        #print("\nPrinting the first 1000 distances inside\n")
        #if len(lsDists) < 1000:
        #    printRng = len(lsDists)
        #else:
        #    printRng = 1000
        #for idx in range(printRng):
        #    print(str(lsDists[idx]))

        totWeight = 0
        sumLongWeightProduct = 0
        # perform the summation of the product of longs with their associated weights
        # as well as the total weight to use in the divisor
        for i in range(len(lsWeights)):
            sumLongWeightProduct += lsLongs[i] * lsWeights[i]
            totWeight += lsWeights[i]
        # last step is to divide the sum of the Long-Weight product by the tot Weight.
        # Check to make sure we don't div by 0
        if totWeight != 0:
            retObj['COM_Long'] = sumLongWeightProduct / totWeight
        # Actually, we need to calculate the area of the county in km^2 so that when
        # calculating the State's COM, we can multiply this county's COM by its appropriate weight.
        # To get this, we multiply our totWeight by the product of the earth's radius squared (in km^2)
        # and the stepAmount (in degrees).
        # Note: we needed to convert from deg to radians in order for an accurate measure. Since
        # this calculation only has one measure in degrees, multiply by one instance of conversion factor.
        retObj['CountyArea'] = totWeight * EarthRad_km**2 * stepAmt * (math.pi/180)
        print("\n\nEstimated area for "+self.name+" county is "+str(retObj['CountyArea'])+" km^2")
        return retObj
    
    def ComputeCOMandArea(self, stepAmt):
        # Returns a dictionary with an ordered pair for the COM and a value for Area
        objOut = {}
        objOut['COM'] = []
        objOut['Area'] = None
        # Need to specify how large step amt will be.
        # We will start with 1 meter increments. At the equator,
        # 1 meter is equivalent to approx 8.983e-6 deg

        # Tells you approx how many distinct lat and long values will be examined
        if stepAmt > 0:
            numLongChecks = int((self.maxLong - self.minLong)/stepAmt)
            numLatChecks = int((self.maxLat - self.minLat)/stepAmt)
            print_s = "\nThis will require approximately "+str(numLongChecks)+" checks along the longitude axis"
            print_s += "\nAnd approximately "+str(numLatChecks)+" checks along the latitude axis"
            print(print_s)
        else:
            print("Must have a step amount > 0")

        # The COM Lat & Long functions return an object containing the respective
        # coordinate as 'COM_Lat' or 'COM_Long' and its area obtained from the computation
        # as 'CountyArea'
        COM_Lat_obj = self.ComputeCOMLat(stepAmt)
        COM_Lat = COM_Lat_obj['COM_Lat']
        CountyArea_fromLat = COM_Lat_obj['CountyArea']
        COM_Long_obj = self.ComputeCOMLong(stepAmt)
        COM_Long = COM_Long_obj['COM_Long']
        CountyArea_fromLong = COM_Long_obj['CountyArea']
        AvgSA = (COM_Lat_obj['CountyArea']+COM_Long_obj['CountyArea']) / 2
        
        objOut['COM'] = [COM_Lat,COM_Long]
        objOut['Area'] = AvgSA

        # set County properties, COM and Area
        self.COM = objOut['COM']
        self.Area = AvgSA
        return objOut

class Segment:
    def __init__(self,coord1,coord2):
        # Attributes helpful for determining whether a coord is within the
        # longitude range and whether it falls to the left, right, or on the segment
        self.coord1 = coord1
        self.coord2 = coord2
        # compute the slope. *Note: if segment is vertical, the slope is
        # undefined, so to denote this, self.slope will be None
        if coord2[1] == coord1[1]:
            self.slope = None
        else:
            self.slope = (coord2[0] - coord1[0])/(coord2[1] - coord1[1])
        # compute intercept (y-intercept, AKA lat-intercept)
        # same caution as with slope. If slope is undefined, intercept = None
        if coord2[1] == coord1[1]:
            self.intercept = None
        else:
            self.intercept = coord1[0] - (self.slope * coord1[1])
    
    def CrossesLong(self,aLong):
        # If the segment is vertical (i.e., Longs are the same b/w coordinates,
        # we do not want to return a Lat because this would
        # invert the set of Lat values computed as within the geographical
        # space at this Longitude.
        if self.coord1[1] == self.coord2[1]:
            return False
        else:
            return (self.coord1[1]<=aLong<=self.coord2[1]) or (self.coord1[1]>=aLong>=self.coord2[1])

    def CrossesLat(self,aLat):
        # If the segment is horizontal (i.e., Lats are the same b/w coordinates,
        # we do not want to return a Lat because this would
        # invert the set of Long values computed as within the geographical
        # space at this Latitude.
        if self.coord1[0] == self.coord2[0]:
            return False
        else:
            return (self.coord1[0]<=aLat<=self.coord2[0]) or (self.coord1[0]>=aLat>=self.coord2[0])


    def getLatAtLong(self,aLong):
        # First check for special cases 1) slope undef, 2) slope is 0
        if self.slope is None:
            # If the segment is vertical (i.e., Longs are the same b/w coordinates),
            # we do not want to return a Lat because this would
            # invert the set of Lat values computed as within the geographical
            # space at this Longitude.
            return None
        elif self.coord1[0]==self.coord2[0]:
            return self.coord1[0]
        else:
            return (self.slope * aLong) + self.intercept
            #return (aLong - self.intercept) / self.slope

    def getLongAtLat(self,aLat):
        # First check for special cases 1) slope undef, 2) slope is 0
        if self.slope == 0:
            # If the segment is horizontal (i.e., slope is 0),
            # we do not want to return a Long because this would
            # invert the set of Long values computed as within the geographical
            # space at this Latitude.
            return None
        elif self.coord1[1] == self.coord2[1]:
            return self.coord1[1]
        else:
            # still double check to make sure we don't cause an error by div by 0
            if self.slope == 0:
                return None
            return (aLat - self.intercept) / self.slope

    def PrintSegmentPoints(self):
        s = str(self.coord1)
        s += '\t-->\t'
        s += str(self.coord2)
        return s

def removeSharedCoords(coords_in):
    arr_coords = coords_in.to_numpy().tolist()
    arr_out_coords = []
    dupl_coords = []
    for row in arr_coords:
        if row in arr_out_coords:
            arr_out_coords.remove(row)
            dupl_coords.append(row)
        elif row not in dupl_coords:
            arr_out_coords.append(row)
    # convert 2D list as a DataFrame object, just as it was provided
    arr_out_coords = pd.DataFrame(arr_out_coords)
    return arr_out_coords

def GetStartingValue(maxVal,stepAmt):
    # want to start at or within the maxLong value, depending
    # on whether the number of decimals in the stepAmt exceeds
    # the maxLong number of decimals or not
    stepAmtDecimals = numDecimals(stepAmt)
    if numDecimals(maxVal) <= stepAmtDecimals:
        startVal = maxVal
    else:
        startVal = round(maxVal,stepAmtDecimals)
        if startVal > maxVal:
            startVal = startVal - pow(10,-(stepAmtDecimals))
    return startVal

def numDecimals(a):
    if type(a) is int:
       return 0
    elif type(a) is float:
        if abs(a) > 1:
            # strip it down to just the decimals (remove value to left of decimal)
            a = a - int(a)
        # take absolute value of a to make positive
        a = abs(a)
        # begin while loop of finding num decimals
        # can only reliably detect up to about 12 decimals
        decCtr = 1
        while abs(a - round(a,decCtr)) > 1e-12 and decCtr <= 11:
            decCtr += 1
        return decCtr
    else:
        return None

def ConvertKMLtoDict(filepath):
    import xml.etree.ElementTree as ET
    import json
    import xmltodict

    with open(filepath) as xml_file:
        data_dict = xmltodict.parse(xml_file.read())

    #for key in data_dict['kml']['Document']:
    #    if key == 'name':
    #        print(data_dict['kml']['Document'][key])
    
    return data_dict

def GetListOfBoundaryCoords(dataIn={}):
    # Receives a dictionary object as parameter from the KML file
    # for the given county.
    # Want to return a list of coordinates, where each coordinate
    # is a 2-number list: [Latitude (in deg), Longitude (in deg)].
    # Note: the first and the last coordinates in the list should
    # be the same, as it has returned to its starting point.
    # 2nd Note: kml files provide data in [longitude,latitude] pairs
    outLs = []
    
    coords = dataIn['kml']['Document']['Placemark']
    for item in coords:
        if 'LineString' in item:
            coordString = item['LineString']['coordinates']
            lsString = coordString.split("\n")
            for idx_lsS in range(len(lsString)):
                lsString[idx_lsS] = lsString[idx_lsS].strip()
                # we should now have a string of just the lat,long values
                # ex., '-82.323978,34.477582'
                if ',' in lsString[idx_lsS]:
                    aCoordStrLs = lsString[idx_lsS].split(',',1)
                    aCoord = [float(strNum) for strNum in aCoordStrLs]
                    # data is provided in opposite order from what we want, so reverse ordered pair
                    aCoord.reverse()
                    outLs.append(aCoord)
    return outLs

def CalculateStateCOM(COMsAreas = {}):
    outObj = {}
    COM = []
    totArea = 0
    WeightedSumLat = 0
    WeightedSumLong = 0
    for c in COMsAreas:
        WeightedSumLat += (c['COM'][0] * c['Area'])
        WeightedSumLong += (c['COM'][1] * c['Area'])
        totArea += c['Area']
    COM.append(WeightedSumLat/totArea)
    COM.append(WeightedSumLong/totArea)
    outObj['COM'] = COM
    outObj['AREA'] = totArea
    return outObj

def ConvertPartialResultsToDict(res_data=""):
    outObj = {}
    # First line contains county name
    # Second line contains COM data: "COM: 32.98812996587019, -81.35775548344009"
    # Third line contains Area data: "SA: 1066.1810320921977"
    linesPerRes = 3
    ls_res = res_data.split('\n')
    numResults = int(len(ls_res)/linesPerRes)
    for i in range(numResults):
        # create object to store COM and AREA
        this_data = {'COM':''
                     ,'Area':''}
        # Get County name and add an object to outObj with county name as key
        resFirstLine = linesPerRes * i
        county_name = str(ls_res[resFirstLine])
        # Make a dict within outObj where key is county_name
        outObj[county_name] = {}

        # Now we need to extract COM and Area data from next two lines
        com_data = str(ls_res[resFirstLine+1])
        sa_data = str(ls_res[resFirstLine+2])
        # Remove the labels 'COM: ' and 'SA: '
        com_data = com_data.replace("COM: ","")
        sa_data = sa_data.replace("SA: ","")
        # Need to find pos of comma in COM data string
        comma_pos = com_data.find(",")
        COM_vals = [float(com_data[:comma_pos]),float(com_data[comma_pos+2:])]
        SA_val = float(sa_data)

        # Write values to county dictionary
        this_data['COM'] = COM_vals
        this_data['Area'] = SA_val
        outObj[county_name] = copy.deepcopy(this_data)

    return outObj

def main():
    ## Code to get list of segments that outline two counties (Berkley and Charleston)
    #csv_read_path = "C:\\Users\\mathl\\OneDrive\Austin's Files\\Personal\\Brainteasers\\SouthCarolinaCentroid\\CountyCSVs\\berkeley_charleston.csv"
    #all_coords = pd.read_csv(csv_read_path,header=None)
    ##print(all_coords.head(10))
    #exterior_coords = removeSharedCoords(all_coords)
    #csv_write_path = "C:\\Users\\mathl\\OneDrive\Austin's Files\\Personal\\Brainteasers\\SouthCarolinaCentroid\\CountyCSVs\\berkeley_charleston_exterior.csv"
    #exterior_coords.to_csv(csv_write_path,header=0,index=0)
    
    ###############################################
    #countyName = "Charleston"
    #c1 = [33.214659, -79.445782]
    #c2 = [33.134197, -79.270001]
    #c3 = [32.493964, -80.330485]
    #c4 = [32.740080, -80.450685]
    #c5 = [33.214659, -79.445782]
    #cList = [c1,c2,c3,c4,c5]
    #aCounty = County(countyName,cList)
    ##print(aCounty.borders[0].PrintSegmentPoints())
    
    #s = aCounty.PrintMyBorders()
    #print(s)

    #aLong = -80
    #lsBorders = aCounty.GetBorderIdxsSpanningLong(aLong)
    #s = "\n\nThe following borders span the given longitude of: "
    #s += str(aLong)
    #s += '\n'+80*'-'
    #for idx in lsBorders:
    #    s += '\nBorder ' + str(idx) + ':\t'
    #    s += aCounty.borders[idx].PrintSegmentPoints()
    #print(s)

    #latObj = aCounty.ComputeCOMLat(1e-5)
    #print(aCounty.name + " county's COM Latitude position is: "+str(latObj['COM_Lat']))
    #print(aCounty.name + " county's total area in km^2 from the COM Latitude method is: "+str(latObj['CountyArea']))
    ###############################################

    # Credits string
    #creditsMatrix = []
    #empty_str = " "
    #line_empty = [empty_str]*80
    #creditsMatrix = [copy.deepcopy(line_empty)]*6
    ##for i in range(80):
    ##    line_empty.append(" ")
    ##for i in range(6):
    ##    creditsMatrix.append(line_empty)
    #print(str(creditsMatrix))
    #for i in range(len(creditsMatrix[0])):
    #    creditsMatrix[0][i] = "~"
    #    creditsMatrix[-1][i] = "~"
    #for i in range(len(creditsMatrix)):
    #    creditsMatrix[i][0] = "~"
    #    creditsMatrix[i][-1] = "~"
    #print(str(creditsMatrix))
    
    text1 = " Code written by Austin T. Hetherington"
    text2 = " All rights reserved. (2022)"

    #i = 0
    #startLoc = 2
    #while i < len(text1):
    #    creditsMatrix[2][i+startLoc] = text1[i]
    #    i += 1
    #print(str(creditsMatrix))
    
    #i = 0
    #startLoc = 2
    #while i < len(text2):
    #    creditsMatrix[3][i+startLoc] = text2[i]
    #    i += 1
    #print(str(creditsMatrix))

    # Convert creditsMatrix to string
    credits = "\n\n"
    #for i in range(len(creditsMatrix)):
    #    for j in range(len(creditsMatrix[i])):
    #        credits += creditsMatrix[i][j]
    #    credits += "\n"
    #print(str(creditsMatrix))
    credits += 80*"~"
    credits += "\n"
    credits += text1 + "\n" + text2
    credits += "\n"
    credits += 80*"~"

    # Create a list to store all Counties
    lsCounties = []
    lsCOMsAreas = []
    stateName = "South Carolina"
    # Set a step amount for the code to increment as it approximates SA
    stepAmt = 1e-5

    # Read the results we've already collected on COM and Area
    #open text file in read mode
    partial_res_filename = "partial_results_extracted_220902.txt"
    text_file = open(partial_res_filename, "r")
 
    #read whole file to a string
    res_data = text_file.read()
 
    #close file
    text_file.close()
 
    partial_res_dict = ConvertPartialResultsToDict(res_data)

    # Loading and parsing a KML file
    folderpath = "C:\\Users\\mathl\\OneDrive\\Austin's Files\\Personal\\Brainteasers\\SouthCarolinaCentroid\\CountyKMLFiles\\All46Counties\\"
    fileNum = 3
    while fileNum <= 48:
        filepath = folderpath + "contourplot (" + str(fileNum) + ").kml"
        kmlDict = ConvertKMLtoDict(filepath)
        # below line gets the county name string from kmlDict keeping only the text before the first comma
        # and then capitalizes the first letter
        countyName = kmlDict['kml']['Document']['name'][:kmlDict['kml']['Document']['name'].find(',')].title()
        countyCoords = GetListOfBoundaryCoords(kmlDict)
        thisCounty = County(countyName, countyCoords)
        #print(thisCounty.PrintMyBorders())
        
        # Check if this county has already been processed; if so, then pull its COM and Area data
        if thisCounty.name in partial_res_dict:
            countyCOMandArea = partial_res_dict[thisCounty.name]
            # Set thisCounty's COM and AREA properties, too
            thisCounty.COM = countyCOMandArea['COM']
            thisCounty.Area = countyCOMandArea['Area']
        else:
            # Otherwise, if it hasn't been computed yet, then we compute
            countyCOMandArea = thisCounty.ComputeCOMandArea(stepAmt)
            # Append these values to the partial_results text file
            # First, make string to append
            appendStr_COMArea = '\n'
            appendStr_COMArea += thisCounty.name
            appendStr_COMArea += '\n'
            appendStr_COMArea += 'COM: '+str(thisCounty.COM[0])+', '+str(thisCounty.COM[1])
            appendStr_COMArea += '\n'
            appendStr_COMArea += 'SA: '+str(thisCounty.Area)
            append_filename = 'partial_results_appending.txt'
            # Open a file with access mode 'a'
            with open(append_filename, "a") as file_object:
                # Append results string at the end of file
                file_object.write(appendStr_COMArea)

        print(thisCounty.StringCOMandArea())
        lsCounties.append(thisCounty)
        lsCOMsAreas.append(countyCOMandArea)
        #print(80*"-")
        fileNum += 1
    
    # Compute wieghted average of all counties COMs by area for answer
    ANSWER = CalculateStateCOM(lsCOMsAreas)
    answerStr = '\n\n'
    answerStr += 80*'*'
    answerStr += '\n\n'
    answerStr += "Based on a step degree amount of "+str(stepAmt)+" we have obtained the following results..."
    answerStr += '\n\n'
    answerStr += 15*'-' + stateName.upper() + 15*'-'
    answerStr += '\n'
    answerStr += "COM:\t"+str(ANSWER['COM'])
    answerStr += '\n'
    answerStr += "AREA:\t"+str(ANSWER['AREA'])+" km^2"
    answerStr += '\n\n'
    answerStr += 80*'*'
    answerStr += '\n\n'

    # Create text file to save computations, answers, etc.
    textFileOutput = ""
    textFileName = stateName.replace(" ","") + "_answer.txt"

    textFileOutput += "Author: Austin T. Hetherington\n"
    textFileOutput += "Date: 2 September 2022\n"
    textFileOutput += "Description: This is the solution to the State center of mass problem.\n"
    textFileOutput += answerStr

    # Show work by printing counties' COM and Area
    for county in lsCounties:
        textFileOutput += county.StringCOMandArea()

    textFileOutput += credits

    try:
        with open(textFileName, 'w') as f:
            f.write(str(textFileOutput))
    except FileNotFoundError:
        print(textFileOutput)

    

main()

