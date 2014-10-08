
import sys
import csv
import traceback
import subprocess as sp
from optparse import OptionParser
import time
import numpy as np

class ScidbQuery:
    """encapsulate the complications calling a SCIDB query: 
     double-quote the query, local or remote server, 
     whether to fetch the result array, the way to output """
    def __init__(self, query=None, server=None, fetch=True, output=None):
        self.query = query
        self.server = server
        self.fetch = fetch
        self.time = 0.0
        self.alreadyRun = False
        self.output = output
        
    def setQuery(self, query):
        "use the same setting to run another query"
        self.query = query
        self.time = 0.0
        self.alreadyRun = False
    
    """ The single quote from SCIDB 13.6 breaks this.  Will reimplement it using pandas when needed   
    def loadToNumpy(self):
        "load the result of a query into a numpy array"
        self.fetch = True
        self.output = sp.PIPE
        sub = self.run( wait=False)
        arr = np.genfromtxt(sub.stdout, delimiter=",", names=True)
        sub.wait()
        return arr
    """
    
    def getRecords(self):
        "return the query as a list of Record"
        self.fetch = True
        self.output = sp.PIPE
        sub = self.run( wait=False)
        results = []
        names = []
        reader = csv.reader(sub.stdout, quotechar="'")
        for row in reader:
            # print row
            if len(names) == 0:
                names = row
            else:
                rec = Record(names)
                rec.setValues(row)
                results.append(rec)
        sub.wait()
        return results
    
    def outputToFile(self, fname=None): 
        "write the output to a file and return the lapsed time"
        self.fetch = True
        if fname is None:
            fname = self.query[0:min(15, len(self.query))]
            fname = fname.translate(None, '()')
        fout = open(fname+'.out', "w")
        self.fetch = True
        self.output = fout
        self.run( wait=True)
        return self.time
    
    def serverActionOnly(self):
        "call to make a change to the server.  No data return is needed.  Return code to indicate success"
        self.fetch = False
        self.output = None
        sub = self.run(wait=True)
        return sub.returncode
    
    def runChain(self, queryChain):
        "run a series of queries with no data returned"
        for i in range(len(queryChain)):
            print 'calling \"%s\" ...' % queryChain[i]
            self.setQuery(queryChain[i])
            if (self.serverActionOnly() != 0 ):
                print '\"%s\" failed' % queryChain[i]
                return -1
        return 0
        
    def run(self, wait=True):
        """the core method to run a query. Take care of all the common complications in one place.
        Most likely you should use one of the convenient methods above"""
        start = time.time()
        if (self.alreadyRun):
            return 0;
        quotedQuery = "'%s'" % self.query
        option = "-aq"
        if not self.fetch:
            option = "-anq"
        if self.server:
            cmd = ['ssh', self.server, '/opt/scidb/13.3/bin/iquery',"-o", "csv+",option, quotedQuery]
        else: # localhost
            cmd = ['iquery',"-o", "csv+", option, self.query]
        sub = sp.Popen(cmd, stdout = self.output, stderr = self.output)
        if wait:
            retcode = sub.wait()
        self.time = time.time() - start
        self.alreadyRun = True
        return sub
           
class Record:
    " A record reconstructed from a csv output file "
    def __init__(self, attrNames):
        self.names = attrNames  
    def setValues(self, attrValues):
        " assume  attrNames and attrValues match up "
        if len(self.names) != len(attrValues):
            print self.names
            print attrValues
            raise Exception("Attribute names and values don't match up.")
        for i in range(len(self.names)):
            setattr(self, self.names[i], attrValues[i])
    def __str__(self):
        result = []
        for name in self.names:
            result.append("%s:%s" % (name, str(getattr(self, name))) )
        return str(result)       
    
def main():
    parser = OptionParser()
    parser.add_option('-q','--query', help='the SCIDB query to call',type="string", dest='query', default='list()')
    (options, args) = parser.parse_args()
    "unit testing"
    sq = ScidbQuery(options.query)
    # sq = ScidbQuery(options.query, server='gtprod31.st-va')
    
    print "testing outputToFile"
    qt = sq.outputToFile()
    print qt, "seconds"
    print "testing getRecords"
    sq.setQuery(options.query)
    recs = sq.getRecords()
    for rec in recs:
        print rec
    print "testing serverActionOnly"
    sq.setQuery(options.query)
    rc = sq.serverActionOnly()
    if rc == 0:
        print "%s succeeded" % options.query
    else:
        print "%s failed" % options.query
        
    """
    print "testing loadToNumpy"
    sq.setQuery(options.query)
    arr = sq.loadToNumpy()
    print arr
    """
    
    sys.exit(); 

if __name__ == "__main__":
    main()
