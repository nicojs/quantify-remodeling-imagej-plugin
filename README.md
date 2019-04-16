# TODO

1. Select 3 images
    * Original file
    * Orientation file
    * Coherency file
2. Provide settings
    * Offset [0,originalImage.width]
    * Border [50 micron]
    * Coherency cutoff [0,1] (.4)
3. Run the script

## How to debug

https://imagej.net/Debugging#Launching_ImageJ_in_debug_mode

1. Compile
    ```
    mvn package
    ```
2. Copy to plugin directory and start ImageJ in debug mode 
    ```
    cp ../../z/github/nicojs/imagejplugins/quantifyremodeling/target/quantifyremodeling-0.0.0-SNAPSHOT.jar plugins/
    ./ImageJ-win64.exe  --debugger=8000
    ```
3. Run your IDE and attach to port 8000
    ```
    -agentlib:jdwp=transport=dt_socket,server=y,suspend=n,address=8000
    ```