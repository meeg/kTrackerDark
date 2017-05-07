#!/usr/bin/env python

import sys
import MySQLdb

sourceFile = sys.argv[1]
sourceSchema = sys.argv[2]
targetSchema = sys.argv[3]

con = MySQLdb.connect(host = 'e906-db1.fnal.gov', port = 3306, user = 'seaguest', passwd = 'qqbar2mu+mu-')
cur = con.cursor()

# get table names
cur.execute("SELECT table_name FROM information_schema.tables WHERE table_schema='%s'" % sourceSchema)
tableNames = [row[0] for row in cur.fetchall()]

# copy the tables first
cur.execute('DROP DATABASE IF EXISTS %s' % targetSchema)
cur.execute('CREATE DATABASE ' + targetSchema)
for table in tableNames:
    cur.execute('CREATE TABLE %s.%s LIKE %s.%s' % (targetSchema, table, sourceSchema, table))
    cur.execute('INSERT %s.%s SELECT * FROM %s.%s' % (targetSchema, table, sourceSchema, table))
con.commit()

# read the new survey
survey = []
for line in open(sourceFile).readlines():
    vals = line.strip().split()
    survey.append(vals)

# make the query
filedNames = survey[0]
for vals in survey[1:]:
	query = 'UPDATE %s.Planes SET ' % targetSchema
	for index, field in enumerate(filedNames):
		if index == 0:
			continue
		query = query + field + '=' + str(vals[index]) + ', '
	query = query[:-2] + " WHERE detectorName='" + vals[0] + "'"

	cur.execute(query)
	print query
con.commit()
