""" io.py - a python module for input/output of objects and data.
-----------------------------------------------------------------------------

(c) Copyright by Paul M. Magwene, 2002-2004  (mailto:pmagwene@sas.upenn.edu)
    Permission to use, copy, modify, and distribute this software and its
    documentation for any purpose and without fee or royalty is hereby granted,
    provided that the above copyright notice appear in all copies and that
    both that copyright notice and this permission notice appear in
    supporting documentation or portions thereof, including modifications,
    that you make.

    THE AUTHOR PAUL M. MAGWENE DISCLAIMS ALL WARRANTIES WITH REGARD TO
    THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
    FITNESS, IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL,
    INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
    FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
    NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
    WITH THE USE OR PERFORMANCE OF THIS SOFTWARE !
-----------------------------------------------------------------------------

Overview

Function            Description
========            ===========
readobject          Read a pickled python object from a file.
writeobject         Pickle a python object ot a file.

tableasrows         Read rows of a delimited file. Return list.
tableascols         Read columns of a delimited file. Returns dict keyed by col indices.
todelimited         Write a table (seq of seqs) to a file.

readseq             Read a endline separated list of strings form a file.
writeseq            Write a seq as a file of endline separate string.

writestr            Write a string to a file.
"""

import os, os.path, cPickle, string, sys
import csv


#------------------------------------------------------------------------------
# filename reformating

def osfilename (fname):
    """Windows style filename string -> unix style.
    """
    fname = string.replace(fname, '\\', '/')
    fname = string.replace(fname, '\t', '/t')
    fname = string.replace(fname, '\n', '/n')
    fname = string.replace(fname, '\r', '/r')
    fname = string.replace(fname, '\f', '/f')
    fname = string.replace(fname, '\v', '/v')
    fname = string.replace(fname, '\a', '/a')
    fname = string.replace(fname, '\b', '/b')
    return os.path.expanduser(fname)

# On platforms other than windows, don't do the filename substitutions    
if 'win' not in sys.platform:
    def osfilename (fname):
        return os.path.expanduser(fname)
        

#------------------------------------------------------------------------------
# Binary I/O based on pickles

def readobject (fname):
    """Routine for loading a pickled (binary) python object."""
    f = open(osfilename(fname),'rb')
    o = cPickle.load(f)
    f.close()
    return o


def writeobject (o, fname):
    """Routine for saving a binary representation of a python object."""
    f = open(os.path.expanduser(fname),'wb')
    cPickle.dump(o,f,1)
    f.close()


#------------------------------------------------------------------------------
# Text I/O in the form of sequences and tables.
                
                    
def tableasrows (filename, delimiter="\t", commentchar="#", skipinitialspace=True,
                hasheader=False, hasrownames=False,
                autoconvert=False, na=False, nastring="NA", navalue=-999):
    """Returns a delimited file as a list of lists.
    
    If hasheader==True, the first row of the file gives the column names, and this
    is returned as a separate list.
    
    If hasrownames==True, the first column of the file give the row names, and this
    is returned as a separate list.
    
    If autoconvert==True, attempts numeric conversions in the order:
    int->float->complex. Otherwise, returns items of table as strings.
    
    If na==True, table cells that contain nastring (default="NA") are replaced
    by the value of navalue (default=-999), and a mask is returned
    indicating the positions of missing values.
    
    Any lines beginning with the commentchar are ignored.
    """
    f = open(osfilename(filename), 'r')
    reader = csv.reader(f, delimiter=delimiter, skipinitialspace=skipinitialspace)
    srows = strip_last_empty(reader, delimiter=delimiter)
    if hasheader:
        header = srows[0]
        srows = srows[1:]
    if hasrownames:
        rownames = []
    rows = []
    mask = []
    for i,row in enumerate(srows):
        try:
            if row[0][0] == commentchar:
                continue
        except IndexError:
            rows.append(row)
            continue
        
        if hasrownames:
            rownames.append(row[0])
            row = row[1:]
                            
        if autoconvert:
            try:    # see if it's an int
                cseq, rowmask = seq_astype(row, type=int, na=na, nastring=nastring, navalue=navalue)
                rows.append(cseq)
                if rowmask is not None:
                    mask.append(rowmask)
                else:
                    mask.append([0]*len(cseq))
            except ValueError:
                try:    # see if it's a float
                    cseq, rowmask = seq_astype(row, type=float, na=na, nastring=nastring, navalue=navalue)
                    rows.append(cseq)
                    if rowmask is not None:
                        mask.append(rowmask)
                    else:
                        mask.append([0]*len(cseq))
                except ValueError:
                    try:    # see if it's a complex number
                        cseq, rowmask = seq_astype(row, type=complex, na=na, nastring=nastring, navalue=navalue)
                        rows.append(cseq)
                        if rowmask is not None:
                            mask.append(rowmask)
                        else:
                            mask.append([0]*len(cseq))
                    except ValueError:  # can't handle it so raise Error
                        print row
                        raise ValueError("Unable to auto-convert table at Row %d. Try again with autoconvert=False" % i)
        else:
            rows.append(row)    

    class TableData(list):
        pass
    
    r = TableData(rows)
    r.data = rows
    r.mask = None
    r.header = None
    r.rownames = None               
    if na:
        r.mask = mask
    if hasheader:
        r.header = header
    if hasrownames:
        r.rownames = rownames
    return r                                  


def tableascols (filename, delimiter="\t", commentchar="#", autoconvert=False):
    """Reads a delimited file and returns cols as as list of lists.
    
    If autoconvert==True, attempts to numeric conversions in order 
    int->float->complex. Otherwise, returns items of table as strings.
    
    Any lines beginning with the commentchar are ignored.
    """
    f = open(osfilename(filename), 'r')
    reader = csv.reader(f, delimiter=delimiter)
    srows = strip_last_empty(reader, delimiter=delimiter)
    columns = {}
    for row in srows:
        if row[0][0] == commentchar:
            continue
        for i, item in enumerate(row):
            if autoconvert:
                try:    # are they ints?
                    columns.setdefault(i, []).append(int(item))
                except ValueError:
                    try:    # are they floats?
                        columns.setdefault(i, []).append(float(item))
                    except ValueError:
                        try:     # are they complex?
                            columns.setdefault(i, []).append(complex(item))
                        except ValueError:  # can't handle conversion
                            raise ValueError("Unable to auto-convert table. Try again with autoconvert=False")
            else:
                columns.setdefault(i, []).append(item)
    f.close()
    
    # now turn dict to list
    keys = columns.keys()
    keys.sort()
    cols = [columns[k] for k in keys]
        
    return cols
    

def seq_astype (seq, type=int, na=False, nastring="NA", navalue=-999):
    """Seq of strings -> Seq w/strings converted to given type.
    """
    try:
        return [type(i) for i in seq], None
    except ValueError:
        if not na:
            raise
        cseq, mask = [],[]
        for i in seq:
            if i == nastring:
                cseq.append(navalue)
                mask.append(1)
            else:
                try:
                    cseq.append(type(i))
                    mask.append(0)
                except ValueError:
                    raise ValueError("Unable to automatically convert seq.")
        return cseq, mask        


def strip_last_empty (rowiter, delimiter="\t"):
    allextra = True
    rows = []
    for row in rowiter:
        if not len(row):
            rows.append(row)
            allextra = False
            continue
        if len(row[-1]):    # if the last element is not an empty string
            allextra = False
        rows.append(row)

    if allextra:            # if every row has the extra empty string, strip it
        return [r[:-1] for r in rows]
    else:
        return rows


def todelimited (seq, filename, header=None, rownames=None, mask=None,
                    delimiter="\t", lineterminator="\n",nastring="NA"):
    """Nested seq -> text file w/elems sep by delimiter, rows sep by endlines.
    """
    f = open(filename, 'w')
    writer = csv.writer(f, delimiter=delimiter, lineterminator=lineterminator)

    newseq = []
    if header is not None:
        writer.writerow(header)
    ct = 0    
    for row in seq:
        nrow = row[:]
        if type(nrow) != type([]):
            nrow = list(nrow)
        if mask is not None:
            for i,each in enumerate(mask[ct]):
                if each: nrow[i] = nastring
        if rownames is not None:
            nrow = [rownames[ct]] + nrow
        writer.writerow(nrow)
        ct += 1
    f.close()
   


def writeseq (v, filename):
    """ sequence -> text file with items of sequence separated by endline.
    """
    f = open(filename, 'w')
    for each in v:
        f.write("%s\n" % str(each))        
    f.close()   


def readseq (filename):
    """ Text file with items of sequence separated by endline -> list of strings.
    """
    f = open(osfilename(filename), 'r')
    items = [string.strip(each) for each in f.xreadlines()]
    f.close()
    return items    



def writestr (s, filename):
    """ Write given string to file.
    """
    f = open(filename, 'w')
    f.write(s)
    f.close()

