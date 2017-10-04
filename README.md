# grafiliv
Grafikkortsliv - life on a graphics card

In order to run a grafiliv simulation:

1.  Obtain binaries by either:

    1.  Compile the source code:

        1.  Download and install the Fluidix library at: <http://www.fluidix.ca/>

        2.  Using the Fluidix app, compile \fluidix\grafiliv\grafiliv.cu

        3.  If you want to compile the GrafilivViewer as well (instead of using the provided binary), download Unity3d <https://unity3d.com/> and open and compile the project at \GrafilivViewer

    2.  Download binaries from the Git repository:

        1.  Download and extract the zip-file corresponding to the version you want to run from <https://github.com/Akodiat/grafiliv/tree/master/app>

2.  Configure \fluidix\grafiliv\config.txt with the parameters you desire.

3.  If compiled for terrain, make sure terrain.stl is present in the same directory as grafiliv.exe

4.  Launch simulation by starting grafiliv.exe. If an earlier simulation was aborted, you have the option to restart it.

5.  Use GrafilivViewer.exe, also located in the same directory as grafiliv.exe, to inspect the simulated organisms