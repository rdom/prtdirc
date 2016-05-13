##Synopsis
```
prtdirc [OPTION] [ARGUMENT] ... [OPTION] [ARGUMENT]

example:
./prtdirc -a 40 -l 0 -x "pi+" -p 1 -w 0 -g 0 -e 1
```
##Options
```

-o    output file name
-i    input file name
-u    look-up file name

-s    run type
                0    simulation
                1    look-up table generation
                2    geometrical reconstruction
                3    likelihood calculation 
                5    calibration
                6    focal plane measurements

-g    geometry configuration
                1    in vacuum
                2    in air
                3    in air + 1cm plastic at front + trigger
                5    for calibration

-gsx    radiator-prism x step in mm

-gsy    radiator-prism y step in mm

-gx  x position of beam-radiator hit point in mm
                0    is a middle of radiator 

-gz  z position of beam-radiator hit point in mm
                0    is at edge of radiator (prism side)

-z    beam dimension
                0    no smearing (default)
                positive int # [mm]   Gaussian smearing in x and y  
                -1  no smearing, random momenta [0,4] GeV


-h    radiator type
                1    bar (default)
                2    plate

-c   MCP layout
                0    MCP covers all FD plain
                2012 
                2014 (default) 

-l    focusing system
                0    no lens
                1    spherical lens
                10   ideal lens (thickness = 0, ideal focusing)

-a    angle between particle beam and bar radiator

-e    number of simulated events
-f     number of the first event for the reconstruction

-x    particle type
              "pi+" 
              "proton"
              "kaon"
                 ...
              "opticalphoton"   

-p    particle momentum [GeV/c]

-w    physical list
                0    standard
                1    without multiple scattering
                10   monochromatic Cherenkov light
                11   10 + 1

-tr    time resolution [ns]
               0.2  (default)  

-r    seed number for the random generator 

-b    batch mode
               1    run silent (without GUI)

```

##Example of script usage from macro folder
```
./ba_scan -j6 -r0 -s5 -e50 -t1 -v0
root da_scan.C'("r_spr39498736070.root","ttt1.root")'
```