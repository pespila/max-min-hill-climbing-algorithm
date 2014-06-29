#!/usr/bin/python

import sys, os
import cgi, cgitb; cgitb.enable()

basedir  = "/var/www/localbene"
workdir  = os.path.join(basedir, "work/work") # cgi must create and delete this
datafile = os.path.join(workdir, "data.txt")
vdfile   = os.path.join(workdir, "data.vd")
dotfile  = os.path.join(workdir, "bene.dot")
netfile  = os.path.join(workdir, "net")
pngfile  = os.path.join(workdir, "bene.png")

bindir   = os.path.join(basedir,"bin")
data2net = os.path.join(bindir,"data2net.sh")
net2parents  = os.path.join(bindir,"net2parents %s -" % netfile)
parents2arcs = os.path.join(bindir,"parents2arcs - -")
arcs2dot     = os.path.join(bindir,"arcs2dot %s - -" % vdfile)
net2dot      = "|".join((net2parents, parents2arcs, arcs2dot))

########### GRAND PLAN               ###########
########### ======================== ###########
########### CREATE WORKDIR OR BUSY   ###########
########### VERIFY INPUT             ###########
########### CREATE VD AND DATA FILES ###########
########### FIND BEST NET            ###########
########### CREATE DOT FILE          ###########
########### CREATE PICTURE           ###########
########### PRINT PICTURE            ###########
########### REMOVE WORKDIR           ###########

def endwork():
    for f in os.listdir(workdir):
        os.remove(os.path.join(workdir, f))
    os.rmdir(workdir)

def perror(title, text, clean=True):
    print """Content-type: text/html

<html>
<body>
<title>Error : %s </title>
<body>
<h1> %s </h1>

<p>
%s
</p>

</body>
</html>
""" % (title, title, text)
    if clean and os.path.exists(workdir):
        endwork()
    sys.exit(0)


########### CREATE WORKDIR OR BUSY   ###########

try:
    os.mkdir(workdir)
except:
    perror("Service Busy !", """
    Somebody else is using the service right now.
    If the problem persists, contact Tomi.
    """, False)
        
########### VERIFY INPUT             ###########

form = cgi.FieldStorage()

essfield = form.getfirst("ess",'').strip()

try:
    if essfield and essfield.upper().endswith("L"):
        ess = float(essfield[:-1])
        essfield = essfield[:-1].strip()+"L"
    else:
        ess = float(essfield)
        
    if ess <= 0:
        perror("Non-positive ESS !",
               """Positive ESS required. (Got "%s".)""" % essfield)

except:
    if not essfield in ("BIC", "AIC", "NML"):
        perror("Illegal ESS-field!",
               """Either positive number (possibly ending with L) or BIC or AIC or NML needed.
               (Got "%s".)""" % essfield)


data_area = form.getfirst("data",'')
datamx = []
for dataline in map(str.strip, data_area.split("\n")):
    if not dataline: continue

    try:
        d = map(int, dataline.split())
    except:
        perror("Error in data format !",
               """The entries on line %d (%s),
               could not be converted to integers. """ % (len(datamx)+1,
                                                          dataline))

    if len(d) > 15:
        perror("Too many variables !",
               """Line %d (%s) has too many (%d) values (i.e. over 15).
               """ % (len(datamx)+1, dataline, len(d)))

    if len(datamx) >= 1000:
        perror("Two many data vectors !",
               """Max 1000 allowed.""")
    
    if min(d) < 0:
        perror("Error in data range !",
               """The entry on line %d (%s)
               is negative. """ % (len(datamx)+1, dataline))
        
    if max(d) > 127:
        perror("Error in data range !",
               """The entry on line %d (%s)
               is greater than 127. """ % (len(datamx)+1, dataline))

    if datamx and len(d) != len(datamx[-1]):
        perror("Error in data format !",
               """Line %d (%s) has different number (%d) of entries 
               than the first  line (that has %d).
               """ % (len(datamx)+1, dataline, len(d), len(datamx[-1])))
    
    datamx.append(d)

if not datamx:
     perror("At least one data vector required !",
               """No data vectors found.
               While not logically impossible, it is no point to try to find
               a best network for zero data vectors.""")
    
    
########### CREATE VD AND DATA FILES ###########

vdf = file(vdfile, "w")
for i, col in enumerate(zip(*datamx)):
    print >>vdf, "V%d\t%s" % (i, "\t".join(map(str, range(max(col)+1))))
vdf.close()

dataf = file(datafile, "w")
for d in datamx:
    print >>dataf, " ".join(map(str, d))                        
dataf.close()

########### FIND BEST NET            ###########

cmd     = " ".join((data2net, vdfile, datafile, essfield, workdir))
cmdin   = os.popen(cmd)
score   = float(cmdin.read())

########### CREATE DOT FILE          ###########

dotf = file(dotfile,"w")
for l in os.popen(net2dot):
    dotf.write(l)
    if l == "  /*@@@BEFORE_NODES@@@*/\n":
        dotf.write('  score [label="Score: %.4f", color=white];\n' % score)
dotf.close()

########### CREATE PICTURE           ###########

os.system("dot -Tpng -o %s %s" % (pngfile, dotfile))

########### PRINT PICTURE            ###########

print "Content-type: image/png\n"
sys.stdout.write(file(pngfile).read())

########### REMOVE WORKDIR           ###########

endwork()
