# Get common recipes
recipes = ENV['fsa'] + '/PlotTools/rake/recipes.rake'
import recipes

require ENV['CMSSW_BASE'] + '/src/FinalStateAnalysis/PlotTools/rake/tools.rb'
require 'pathname'

$jobid = ENV['jobid']
$blind = ENV['blind']

# Figure out what run period we are in
$period = '13TeV'
PU = ENV['PU']
#if $jobid.include? '8TeV'
#  $period = '8TeV'
#end


################################################################################
## Sample names ################################################################
################################################################################
#
# Get sample names containing a substring
def get_sample_names(substring)
  inputs = Dir.glob("inputs/#{$jobid}/*.txt").select {|x| x.include? substring}
  inputs = inputs.map{|x| File.basename(x).sub(".txt", "")}
  return inputs
end
#

samples = Hash[
                "run2017Av2" => get_sample_names('data_SingleMuon_Run2017A_PromptReco-v2'),
                "run2017Av3" => get_sample_names('data_SingleMuon_Run2017A_PromptReco-v3'),
                "run2017Bv1" => get_sample_names('data_SingleMuon_Run2017B_PromptReco-v1'),
     		"run2017Bv2" => get_sample_names('data_SingleMuon_Run2017B_PromptReco-v2'),
		"run2017Cv1" => get_sample_names('data_SingleMuon_Run2017C_PromptReco-v1'),
                "run2017Cv2" => get_sample_names('data_SingleMuon_Run2017C_PromptReco-v2'),
                "run2017Cv3" => get_sample_names('data_SingleMuon_Run2017C_PromptReco-v3'),
     		"run2017Dv1" => get_sample_names('data_SingleMuon_Run2017D_PromptReco-v1'),
		"run2017Ev1" => get_sample_names('data_SingleMuon_Run2017E_PromptReco-v1'),
		"wmHToTT" => get_sample_names('WminusHToTauTau_M125_13TeV_powheg_pythia8_v6'),
        	"wpHToTT" => get_sample_names('WplusHToTauTau_M125_13TeV_powheg_pythia8_v6'),
        	"wgLNuEE" => get_sample_names('WGstarToLNuEE_012Jets_13TeV-madgraph_v6-v1'),
        	"wgLNuMuMu" => get_sample_names('WGstarToLNuMuMu_012Jets_13TeV-madgraph_v6-v1'),
        	"zzTo2L2Q" => get_sample_names('ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1'),
        	"zzTo4L" => get_sample_names('ZZTo4L_13TeV-amcatnloFXFX-pythia8_v6_ext1-v1'),
        	"ewkwm2JWToLNu" => get_sample_names('EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6'),
        	"ewkwp2JWToLNu" => get_sample_names('EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6'),
        	"ewkz2JZToLL" => get_sample_names('EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v6'),
        	"ewkz2JZToNuNu" => get_sample_names('EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v6'),
        	"wGToLNuG" => get_sample_names('WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6'),
        	"ww" => get_sample_names('WW_TuneCUETP8M1_13TeV-pythia8_v6'),
        	"wz" => get_sample_names('WZ_TuneCUETP8M1_13TeV-pythia8_v6'),
        	"zz" => get_sample_names('ZZ_TuneCUETP8M1_13TeV-pythia8_v6'),
        	"wzJToLLLNu" => get_sample_names('WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8_v6-v1'),
        	"wzTo1L1Nu2Q" => get_sample_names('WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v3'),
        	"wzTo1L3Nu" => get_sample_names('WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1'),
        	"wzTo2L2Q" => get_sample_names('WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1'),
        	"wwTo1L1Nu2Q" => get_sample_names('WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1'),
        	"www4F" => get_sample_names('WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_v6-v1'),
        	"dy1JToLLM10to50" => get_sample_names('DY1JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1'),
        	"dy1JToLLM50" => get_sample_names('DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1'),
        	"dy2JToLLM10to50" => get_sample_names('DY2JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1'),
        	"dy2JToLLM50" => get_sample_names('DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1'),
        	"dy3JToLLM50" => get_sample_names('DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1'),
        	"dy4JToLLM50" => get_sample_names('DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1'),
        	"dyJToLLM10to50" => get_sample_names('DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1'),
        	"dyJToLLM50" => get_sample_names('DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6'),
        	"w1JetsToLNu" => get_sample_names('W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1'),
        	"w2JetsToLNu" => get_sample_names('W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6'),
        	"w3JetsToLNu" => get_sample_names('W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6'),
        	"w4JetsToLNu" => get_sample_names('W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6'),
        	"wJetsToLNu" => get_sample_names('WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6'),
		"ttH" => get_sample_names('ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8_v6_ext4-v1'),
		"TT" => get_sample_names('TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_v6-v1'),
		"TTevtgen" => get_sample_names('TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen_v6-v1'),
		"STantitop4f" => get_sample_names('ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_v6-v1'),
		"STtop4f" => get_sample_names('ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_v6-v1'),
		"STtWantitop5f" => get_sample_names('ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v6_ext1-v1'),
		"STtWtop5f" => get_sample_names('ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v6_ext1-v1')
]



# Function to get the .root files for an analyzer and samples
def get_analyzer_results(analyzer, the_samples)
  output = Array.new
  analyzer_base = analyzer.sub('.py', '')
  the_samples.each do |sample|
    output << "results/#{$jobid}/#{analyzer_base}/#{sample}.root"
  end
  return output
end
def get_plotter_results(analyzer)
  output = Array.new
  analyzer_base = analyzer.sub('.py', '')
end

#task :MuMu => get_analyzer_results("MuMuAnalyser.py", samples['run2017Bv1'] + samples['run2017Bv2'] + samples['run2017Cv1'] + samples['run2017Cv2'] + samples['run2017Cv3'] + samples['run2017Dv1'])

#task :MuMu => get_analyzer_results("MuMuAnalyser.py", samples['wmHToTT'] + samples['wpHToTT'] + samples['wgLNuEE'] + samples['wgLNuMuMu'] + samples['zzTo2L2Q'] + samples['zzTo4L'] + samples['ewkwm2JWToLNu'] + samples['ewkwp2JWToLNu'] + samples['ewkz2JZToLL'] + samples['ewkz2JZToNuNu'] + samples['wGToLNuG'] + samples['ww'] + samples['wz'] + samples['zz'] + samples['wzJToLLLNu'] + samples['wzTo1L1Nu2Q'] + samples['wzTo1L3Nu'] + samples['wzTo2L2Q'] + samples['wwTo1L1Nu2Q'] + samples['www4F'] + samples['dy1JToLLM10to50'] + samples['dy1JToLLM50'] + samples['dy2JToLLM10to50'] + samples['dy2JToLLM50'] + samples['dy3JToLLM50'] + samples['dy4JToLLM50'] + samples['dyJToLLM10to50'] + samples['dyJToLLM50'] + samples['w1JetsToLNu'] + samples['w2JetsToLNu'] + samples['w3JetsToLNu'] + samples['w4JetsToLNu'] + samples['wJetsToLNu'] + samples['ttH'] + samples['TT'] + samples['TTevtgen'] + samples['STantitop4f'] + samples['STtop4f'] + samples['STtWantitop5f'] + samples['STtWtop5f'])

#task :MM => get_analyzer_results("MuMuAnalyser.py", samples['run2017Bv1'] + samples['run2017Bv2'] + samples['run2017Cv1'] + samples['run2017Cv2'] + samples['run2017Cv3'] + samples['run2017Dv1'] + samples['wmHToTT'] + samples['wpHToTT'] + samples['wgLNuEE'] + samples['wgLNuMuMu'] + samples['zzTo2L2Q'] + samples['zzTo4L'] + samples['ewkwm2JWToLNu'] + samples['ewkwp2JWToLNu'] + samples['ewkz2JZToLL'] + samples['ewkz2JZToNuNu'] + samples['wGToLNuG'] + samples['ww'] + samples['wz'] + samples['zz'] + samples['wzJToLLLNu'] + samples['wzTo1L1Nu2Q'] + samples['wzTo1L3Nu'] + samples['wzTo2L2Q'] + samples['wwTo1L1Nu2Q'] + samples['www4F'] + samples['dyJToLLM50'] + samples['wJetsToLNu'] + samples['ttH'] + samples['TT'] + samples['TTevtgen'] + samples['STantitop4f'] + samples['STtop4f'] + samples['STtWantitop5f'] + samples['STtWtop5f'])

task :MM => get_analyzer_results("MuMuAnalyser.py", samples['dy1JToLLM10to50'] + samples['dy2JToLLM10to50'] + samples['dyJToLLM10to50'] + samples['dy1JToLLM50'] + samples['dy2JToLLM50'] + samples['dy3JToLLM50'] + samples['dy4JToLLM50'] + samples['w1JetsToLNu'] + samples['w2JetsToLNu'] + samples['w3JetsToLNu'] + samples['w4JetsToLNu'])

#samples['run2017Av2'] + samples['run2017Av3'] + samples['run2017Bv1'] + samples['run2017Bv2'] + samples['run2017Cv1'] + samples['run2017Cv2'] + samples['run2017Cv3'] + samples['run2017Dv1'] + samples['run2017Ev1']


$etdir = "plots/#{$jobid}/MuMuAnalyzer/mt_now/"
directory $etdir 
file  "#{$etdir}/plot#{$period}.root" do |t|
  sh "echo $jobid"
  sh "python myNewPlotterReco.py" 
  
end