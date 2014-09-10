# Loads chromosome differences from Wormbase

import glob
from peewee import *


db = SqliteDatabase('chrom_diff.db')
db.connect()

# Define Classes
class BaseModel(Model):
    class Meta:
        database = db

class liftover(BaseModel):
    build = IntegerField()
    organism = CharField()
    chromosome = CharField()
    mismatch_start1 = IntegerField()
    mismatch_end1 = IntegerField()
    length1 = IntegerField()
    mismatch_start2 = IntegerField()
    mismatch_end2 = IntegerField()
    length2 = IntegerField()
    flipped = BooleanField()


#=====================#
# Load Wormbase Files #
#=====================#

# Remove tables and reinitialize
db.drop_tables([liftover], safe=True)
db.create_tables([liftover], safe=True)

wb_files = glob.glob("Chrom_Diff/sequence_diff*")

for f in wb_files:
	print f
	# Get Wormbase ID
	wormbase_id = f.title().replace("Chrom_Diff/Sequence_Differences.Ws","")
	# Parse File
	for line in open(f, 'r').readlines():
		if line.startswith("Chromosome"):
			cur_chrom = line.split(" ")[1].replace("\n","")
		else:
			# Skip blank lines
			if line.startswith(" ") == False:
				coordinates = [int(x) for x in line.replace("\n","").split("\t")]
				lf = liftover()
				lf.build = wormbase_id
				lf.organism = "C. elegans"
				lf.chromosome = cur_chrom
				# Add 1 for 1-based coordinates.
				lf.mismatch_start1 = coordinates[0] + 1 
				lf.mismatch_end1 = coordinates[1] + 1
				lf.length1 = coordinates[2]
				lf.mismatch_start2 = coordinates[3] + 1
				lf.mismatch_end2 = coordinates[4] + 1
				lf.length2 = coordinates[5]
				if len(line) == 6:
					lf.flipped = coordinates[6]
				else:
					lf.flipped = 0
				lf.save()
