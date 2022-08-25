#!/usr/bin/ruby
# == NAME
# covid_pipeline.rb
#
# == USAGE
# ./this_script.rb [ -h | --help ]
#[ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# A skeleton script for Ruby
#
# == OPTIONS
# -h,--help Show help
# -i,--infile=INFILE input file
# -o,--outfile=OUTFILE : output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#  Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'
require 'open3'

### Define modules and classes here

def msg(string)

	warn "#{Time.now}: #{string.strip}"

end
def cmd(string)

	msg("Executing: #{string}")
	system(string)

end

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "A script description here"
opts.separator ""
opts.on("-f","--folder", "=FOLDER","Input folder") {|argument| options.folder = argument }
opts.on("-x","--force", "Force creating of output folder") {|argument| options.force = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 


###################################
INPUT_FOLDER = options.folder
INPUT_CELL = options.folder.split("/")[-1]
RUN_DATE = options.folder.split("/")[-1].split("_")[0]
SAMPLE_DATE = RUN_DATE.to_i-2
year,month,day = SAMPLE_DATE.to_s.scan(/.{1,2}/)
SAMPLE_DATE_FOLDER = "20#{year}_#{month}_#{day}"
#OUTPUT_BASE = "/mnt/diagx_out/sarscov2"
OUTPUT_BASE = "/work_ifs/sukmb352/projects/covid/pipeline_test"
OUTDIR = "#{OUTPUT_BASE}/#{SAMPLE_DATE_FOLDER}/CorSurV_Run_#{INPUT_CELL}"
###################################

msg("--------------------------------------------")
msg("Welcome to the Sars-Cov2 Pipeline Automator, dear #{ENV['USER']}!")
msg("--------------------------------------------")
msg("")
msg("Compiling analysis environment...")
msg("Input folder: 	#{options.folder}")
msg("Run date: 		#{RUN_DATE}")
msg("Sample date:		#{SAMPLE_DATE_FOLDER}")
msg("Result dir:		#{OUTDIR}")
msg("--------------------------------------------")

# Check input folder

raise "Specified input folder not found (#{options.folder})" unless Dir.exist?(options.folder)

msg("Checking sanity of the input data...")
## Validation ###

expected_folders = [ "DX_CorSurV_WGS", "ngs_QC-Samples_01" ]

expected_folders.each do |ef|
	raise "Did not find expected sub-folder #{ef} in the target folder" unless Dir.exist?("#{options.folder}/#{ef}")
end

#################

msg("...Checking target directory...")
raise "Missing the target directory!" unless Dir.exist?(OUTPUT_BASE)
msg("...Checking planned result directory...")
raise "Result directory already exists - use -x to overwrite it" if Dir.exist?(OUTDIR) && !options.force
msg("...Creating the result directories...")
cmd("mkdir -p #{OUTDIR}")
expected_folders.each do |ef|
	cmd("mkdir -p #{OUTDIR}/#{ef}")
end
msg("Updating pipeline code...")
Open3.popen3("nextflow pull ikmb/virus-pipe") do |stdin,stdout,stderr,wait_thr|
	pid = wait_thr.pid
	exit_status = wait_thr.value
	msg("Pipeline pull reporting: #{exit_status}")
end
msg("Starting pipeline run(s) - this will take a while...")
expected_folders.each do |ef|
	msg("...#{ef}")

	input = "#{INPUT_FOLDER}/#{ef}"
	results = "#{OUTDIR}/#{ef}/results"
	
	nfx_cmd = "nextflow run ikmb/virus-pipe --folder #{input} --outdir #{results} --run_name #{SAMPLE_DATE} -profile diagnostics"
	cmd(nfx_cmd)
	cmd("rm -Rf work")
end

msg("Pipeline runs completed")

msg("Time to make a pretty report....!")
Dir.chdir(OUTPUT_BASE) do |dir|
	cmd("ruby report.rb")
	msg("Report generated in #{OUTPUT_BASE}")
end
