
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class getmissingjobs {
    String parameters[] = new String[3];

    public getmissingjobs (int minDirectory, int maxDirectory, int minPoint, int maxPoint, int maxJobID, String N, String gaps, 
            String G, String mN, String mS, String mD, String startArgumentsFileName, String fileRunningJobsName, String filterText, 
            String fileJobIDName, int startPreviousIteration, String[] fileDirectoriesName, String[] filePointsName, String[] fileLambdaIndexesName, int minLambdaIndex, int maxLambdaIndex) {
        try {
            parameters[0]=fileDirectoriesName[1];
            parameters[1]=filePointsName[1];
            parameters[2]=fileLambdaIndexesName[1];
            
            int licznikStartArguments=0;
            BufferedReader startArguments = new BufferedReader(new FileReader(startArgumentsFileName));
            String buffer; while ((buffer=startArguments.readLine())!=null) licznikStartArguments++; startArguments.close(); 
            String[] startArgumentsList = new String[licznikStartArguments]; 
            startArguments = new BufferedReader(new FileReader(startArgumentsFileName));
            for (int i=0;i<licznikStartArguments;i++) {
                buffer = startArguments.readLine();
                startArgumentsList[i] = buffer.split("\t")[0];
            }
            startArguments.close();
            
            boolean[][][] jobFailedList = new boolean[maxLambdaIndex][maxDirectory][maxPoint];
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
                           saveLambdaIndexes = new BufferedWriter(new FileWriter(fileLambdaIndexesName[0],false)),
                           saveJobIDs = new BufferedWriter(new FileWriter(fileJobIDName,false));
            for (int l=minLambdaIndex;l<jobFailedList.length;l++) for (int i=minDirectory-1;i<jobFailedList[0].length;i++) for (int j=minPoint;j<jobFailedList[0][0].length;j++) if (jobFailedList[l][i][j]) {
                saveDirectories.write(String.valueOf(i+1)); saveDirectories.newLine();
                savePoints.write(String.valueOf(j)); savePoints.newLine();
                saveLambdaIndexes.write(String.valueOf(l)); saveLambdaIndexes.newLine();

                String actualFolderList[] = new File("2D_N-"+N+"_gaps-"+gaps+"_G-"+G+"_badanie-"+String.valueOf(i+minDirectory)+"_mN-"+mN+"_mS-"+mS+"_mD-"+mD).list();
                int minJobID = maxJobID;
                if (actualFolderList!=null) for (int k=0;k<actualFolderList.length;k++) 
                    if (actualFolderList[k].contains(startArgumentsList[G.equals("1")?j:(licznikStartArguments-1-j)]) && actualFolderList[k].contains("_l-"+String.valueOf(l))) 
                        if (actualFolderList[k].contains("Configurations") && !actualFolderList[k].endsWith("Results.txt")) 
                            try{minJobID = Math.min(minJobID,Integer.parseInt(actualFolderList[k].substring(2).split("_")[0]));}catch(Exception e1){}
                saveJobIDs.write(String.valueOf(minJobID+startPreviousIteration)); saveJobIDs.newLine();
            }
            saveDirectories.close(); savePoints.close(); saveLambdaIndexes.close(); saveJobIDs.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void main(String[] args) {
        new getmissingjobs(Integer.parseInt(args[0]),Integer.parseInt(args[1]),Integer.parseInt(args[2]),Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5],args[6], 
            args[7],args[8],args[9],args[10],args[11],args[12],args[13],args[14],Integer.parseInt(args[15]),args[16].split("!"),args[17].split("!"),args[18].split("!"),Integer.parseInt(args[19]),Integer.parseInt(args[20]));
    }
    
}