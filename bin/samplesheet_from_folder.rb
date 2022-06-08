#!/bin/env ruby

require 'optparse'
require 'ostruct'

### Define modules and classes here

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "Reads Fastq files from a folder and writes a sample sheet to STDOUT"
opts.separator ""
opts.on("-f","--folder", "=FOLDER","Folder to scan") {|argument| options.folder = argument }
opts.on("-c","--centre", "=CENTRE","Name of sequencing centre") {|argument| options.centre = argument }
opts.on("-p","--platform", "=PLATFORM","Name of sequencing instrument") {|argument| options.platform = argument }
opts.on("-s","--sanity", "Perform sanity check of md5 sums") { options.sanity = true }
opts.on("-l","--lookup", "=LOOKUP","Add patient IDs from lookup table") { |argument| options.lookup = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 

abort "Folder not found (#{options.folder})" unless File.directory?(options.folder)

date = Time.now.strftime("%Y-%m-%d")
options.centre ? center = options.centre : center = "IKMB"

fastq_files = Dir["#{options.folder}/*_R*.fastq.gz"]

groups = fastq_files.group_by{|f| f.split("/")[-1].split(/_S/)[0] }

warn "Building input sample sheet from FASTQ folder"
warn "Performing sanity check on md5sums" if options.sanity

lookup_table = {}
if options.lookup
	warn "Using lookup table to attach patient IDs"
	lines = IO.readlines(options.lookup)
	lines.each do |line|
		next unless line.include?(",")
		lib,patient = line.strip.split(",")
		lookup_table[lib] = patient
	end
end
	
options.platform ? sequencer = options.platform : sequencer = "NovaSeq6000"

puts "IndivID;SampleID;RGID;R1;R2"

patients_replaced = 0
# group = the library id, may be split across lanes
groups.each do |group, files|

	warn "...processing library #{group}"

	pairs = files.group_by{|f| f.split("/")[-1].split(/_R[1,2]/)[0] }

	pairs.each do |p,reads|

        	left,right = reads.sort.collect{|f| File.absolute_path(f)}
	
		abort "This sample seems to not be a set of PE files! #{p}" unless left && right

		if options.sanity
			Dir.chdir(options.folder) {
				[left,right].each do |fastq|
					fastq_simple = fastq.split("/")[-1].strip
					raise "Aborting - no md5sum found for fastq file #{fastq}" unless File.exists?(fastq_simple + ".md5")
					status = `md5sum -c #{fastq_simple}.md5`
					raise "Aborting - failed md5sum check for #{fastq}" unless status.strip.include?("OK")
				end
			}
		end

		# H26247-L3_S1_L001_R1_001_fastqc.html
        	library = group.split("_S")[0]
        	sample = group.split("_S")[0]
		if group.match(/^S.*_K[0-9]*.*/)
			sample = group.split("_")[1..-1].join("_").split("_S")[0]
		end

		individual = group.split("-")[0]

		if lookup_table.has_key?(sample)
			individual = lookup_table[sample]
			patients_replaced += 1
		end

        	e = `zcat #{left} | head -n1 `
		header = e

        	instrument,run_id,flowcell_id,lane,tile,x,y = header.split(" ")[0].split(":")

		index = header.split(" ")[-1].split(":")[-1]
        	readgroup = flowcell_id + "." + lane + "." + library 

        	pgu = flowcell_id + "." + lane + "." + index

        	puts "#{individual};#{sample};#{readgroup};#{left};#{right}"
	end
end

if (options.lookup)
	warn "Replaced #{patients_replaced} patient IDs from lookup file"
end
