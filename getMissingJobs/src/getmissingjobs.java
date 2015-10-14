
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class getmissingjobs {
    String parameters[] = new String[3];

    public getmissingjobs (int minDirectory, int maxDirectory, int minPoint, int maxPoint, int maxJobID, String N, String gaps, 
            String G, String mN, String mD, String startArgumentsFileName, String fileRunningJobsName, String filterText, 
            String fileJobIDName, int startPreviousIteration, String[] fileDirectoriesName, String[] filePointsName, String[] filePressureIndexesName) {
        try {
            parameters[0]=fileDirectoriesName[1];
            parameters[1]=filePointsName[1];
            parameters[2]=filePressureIndexesName[1];
            
            int licznik=0;
            BufferedReader startArguments = new BufferedReader(new FileReader(startArgumentsFileName));
            String buffer; while ((buffer=startArguments.readLine())!=null) licznik++; startArguments.close(); 
            String[][] startArgumentsList = new String[licznik][2]; 
            startArguments = new BufferedReader(new FileReader(startArgumentsFileName));
            for (int i=0;i<licznik;i++) {
                buffer = startArguments.readLine();
                for (int j=0;j<2;j++) startArgumentsList[i][j] = buffer.split("\t")[j];
            }
            startArguments.close();
            
            boolean[][][] jobFailedList = new boolean[licznik][maxDirectory][maxPoint];
            for (int i=0;i<jobFailedList.length;i++) for (int j=0;j<jobFailedList[0].length;j++) for (int k=0;k<jobFailedList[0][0].length;k++) jobFailedList[i][j][k]=true;
            
            BufferedReader runningJobs = new BufferedReader(new FileReader(fileRunningJobsName));
            while ((buffer=runningJobs.readLine())!=null) {
                String[] result = buffer.split(" ");
                for (int i=0;i<result.length;i++) if (result[i].contains(filterText)) {
                    result[i] = result[i].substring(filterText.length());
                    int parametersValues[] = new int[parameters.length];
                    for (int j=0;j<parameters.length;j++) {
                        int parameterStartIndex = result[i].indexOf(parameters[j]), parameterEndIndex=parameterStartIndex+1;
                        if (parameterStartIndex>=0) {
                            while (parameterEndIndex<result[i].length() && Character.isDigit(result[i].charAt(parameterEndIndex))) parameterEndIndex++;
                            parametersValues[j]=Integer.parseInt(result[i].substring(parameterStartIndex+1,parameterEndIndex));
                        }
                    }
                    try {jobFailedList[parametersValues[2]][parametersValues[0]-1][parametersValues[1]]=false;} catch (Exception e1) {}  //moga byc parametry WYZSZE [outOfBounds] zadane w innej serii
                    break;
                }
            }
            runningJobs.close();
            
            BufferedWriter saveDirectories = new BufferedWriter(new FileWriter(fileDirectoriesName[0],false)),
                           savePoints = new BufferedWriter(new FileWriter(filePointsName[0],false)),
                           savePressureIndexes = new BufferedWriter(new FileWriter(filePressureIndexesName[0],false)),
                           saveJobIDs = new BufferedWriter(new FileWriter(fileJobIDName,false));
            for (int l=0;l<licznik;l++) for (int i=minDirectory-1;i<jobFailedList[0].length;i++) for (int j=minPoint;j<jobFailedList[0][0].length;j++) if (jobFailedList[l][i][j]) {
                saveDirectories.write(String.valueOf(i+1)); saveDirectories.newLine();
                savePoints.write(String.valueOf(j)); savePoints.newLine();
                savePressureIndexes.write(String.valueOf(l)); savePressureIndexes.newLine();

                String actualFolderList[] = new File("2D_N-"+N+"_gaps-"+gaps+"_G-"+G+"_badanie-"+String.valueOf(i+1)+"_mN-"+mN+"_mD-"+mD+"_p-"+startArgumentsList[G.equals("1")?l:(licznik-1-l)][1]).list();
                int minJobID = maxJobID;
                if (actualFolderList!=null) for (int k=0;k<actualFolderList.length;k++) 
                    if (actualFolderList[k].contains(startArgumentsList[G.equals("1")?j:(licznik-1-j)][0])) 
                        if (actualFolderList[k].contains("Configurations") && !actualFolderList[k].endsWith("Results.txt")) 
                            try{minJobID = Math.min(minJobID,Integer.parseInt(actualFolderList[k].substring(2).split("_")[0]));}catch(Exception e1){}
                saveJobIDs.write(String.valueOf(minJobID+startPreviousIteration)); saveJobIDs.newLine();
            }
            saveDirectories.close(); savePoints.close(); savePressureIndexes.close(); saveJobIDs.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void main(String[] args) {
        new getmissingjobs(Integer.parseInt(args[0]),Integer.parseInt(args[1]),Integer.parseInt(args[2]),Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5],args[6], 
            args[7],args[8],args[9],args[10],args[11],args[12],args[13],Integer.parseInt(args[14]),args[15].split("!"),args[16].split("!"),args[17].split("!"));
    }
    
}