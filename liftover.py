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



# Remap Coordinates
def shift_pos(pos, pos_type, mis_start1, mis_end1, length1, mis_start2, mis_end2, length2, shift_direction):
	# Shift Right if beyond end position
	if ((length1 > 0 and pos >= mis_end1) or
		(length1 == 0 and pos > mis_end1)):
		return pos + ((length2 - length1) * shift_direction)
	elif (pos >= mis_start1):
	    # within a changed segment; if the position within the original segment
	    # maps to beyond the end of the replacement segment, then this position no
	    # longer exists, in which case we barf
	    #print pos - mis_start1, length2
	    if (pos - mis_start1 + 1 > length2): 
	    	print "INVOKED"
	    	#return pos + ((length2 - length1) * shift_direction)
	    	return (mis_start2 + (length2*shift_direction))
	    else:
	    	return pos
	else:
		return pos

def remap_coords(organism, release1, release2, chromosome, start, end = None):

	# Determine direction
	if (release1 > release2):
		# Reverse liftover
		liftover_direction = liftover.build.desc()
		release1, release2 = release2, release1
		shift_direction = -1
	else:
		# Forward liftover
		liftover_direction = liftover.build.asc()
		shift_direction = 1

	# Read Mapping Data
	mapping_coords = liftover.select().where(liftover.build >= release1, liftover.build <= release2, liftover.chromosome==chromosome).order_by(liftover_direction).dicts()
	for i in mapping_coords:
		
		cur_length = end - start

		start = shift_pos(start, "start", i["mismatch_start1"], i["mismatch_end1"], i["length1"], i["mismatch_start2"], i["mismatch_end2"],  i["length2"], shift_direction)
		end = shift_pos(end, "end", i["mismatch_start1"], i["mismatch_end1"], i["length1"], i["mismatch_start2"], i["mismatch_end2"],  i["length2"], shift_direction)


		indel = (cur_length,(end-start))

		#print start, end, start - end, i["build"]
	return start, end, indel

def report_error(mapped, orig, pos_type, direction, k, v):
	if (int(mapped) != orig):
		print "%2s Remap %5s Fail - %10s %10s - %10s - %s" % (k+1, pos_type, mapped, orig, v["old_chrom"], direction) 
		return 1
	else:
		return 0


def debug():
	import csv
	test_gff = csv.DictReader(open("test/test_coords2.ws75.ws76.txt"), delimiter="\t")
	fail_count = 0
	for k,v in enumerate(test_gff):
		#print k,v[0], v[1], v[5]

		# Forward
		start, end, indel = remap_coords("C. elegans", 75, 76, v["old_chrom"], int(v["old_start"]), int(v["old_end"]))
		fail_count += report_error(v["new_start"],start, "start", "Forward", k,v)
		fail_count += report_error(v["new_end"],end, "end", "Forward", k,v)

		# Backwards
		start, end, indel = remap_coords("C. elegans", 76, 75, v["old_chrom"], int(v["new_start"]), int(v["new_end"]))
		fail_count += report_error(v["old_start"],start, "start", "Reverse", k, v)
		fail_count += report_error(v["old_end"],end, "end", "Reverse", k, v)

	print "Total Failures: %s / %s" % (fail_count, k)

			#print "End remap failed: %s - %s" % (v["old_end"], end)

debug()

remap_coords("C. elegans", 220, 235,"CHROMOSOME_IV", 9403845, 9403855)
remap_coords("C. elegans", 220, 235,"CHROMOSOME_IV", 9403845, 9403855)