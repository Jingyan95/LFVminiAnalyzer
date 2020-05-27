# Top LFV mini Analyzer  

This is a package consisting of one analyzer and three root files (MC samples) produced by Reza with https://github.com/rgoldouz/miniAOD_IIHE_2018

## I. File Lists

<table border="0">
 <tr>
    <td><b style="font-size:30px">Files</b></td>
    <td><b style="font-size:30px">Description</b></td>
 </tr>
 <tr>
    <td>LFVAnalyzer.C</td>
    <td>A .C file we use to analyze samples</td>
 </tr>
 <tr>
    <td>outfile_TTTo2L2Nu.root</td>
    <td>Simulated background MC sample</td>
 </tr>
 <tr>
    <td>outfile_WZTo3LNu.root</td>
    <td>Simulated background MC sample</td>
 </tr>
 <tr>
    <td>STJets_13TeV_LFV_utemu_vector_Madgraph_outfile.root</td>
    <td>Simulated Signal MC sample</td>
 </tr>
</table>

## II. To open a root file

```sh
git clone https://github.com/Jingyan95/LFVminiAnalyzer.git
cd LFVminiAnalyzer
root -l STJets_13TeV_LFV_utemu_vector_Madgraph_outfile.root
```

This opens the signal sample as a pointer _file0, then you can do:

```sh
_file0->ls()
This lists the contents of the root file. 
```

In this ntuple, there are two different objects: TTree and TH1F, TH1F is simply histogram which you can draw with:

```sh
pileupDist->Draw()
This is probably slow if you are in the US.
```

TTree is a collection of TBranch objects which are containers we use to store an array of data. You can do:

```sh
IIHEAnalysis->Print()
This will print all the branches that are stored in the TTree.

IIHEAnalysis->Draw("gsf_pt")
This draws the distribution of the gsf electron pt.

IIHEAnalysis->Scan("gsf_pt")
This prints out some data. The "scan" command should be fast. 
```

## III. To run/modify LFVAnalyzer.C

```sh
root -l
.x LFVAnalyzer.C
This produces two histograms with the extension pdf
```

There are many ways to modify LFVAnalyzer.C, for example, we can use vim editor:

```sh
vim LFVAnalyzer.C
enter "I" to turn on "edit" mode and enter "ese" to go back to "readonly" mode.
After you are done with editing, enter ":w" to save the changes and ":q" to exit.
```


