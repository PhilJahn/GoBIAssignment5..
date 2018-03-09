import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class GOEnrich {

	public static void main(String[] args) {
		
		String countsfilesPath ="";
		String labelsPath ="";
		String outputPath ="";
		String configPath = "";
		for(int i =0; i < args.length-1; i++){
			if(args[i].equals("-countfiles")){
				countsfilesPath = args[i+1];
				i++;
			}
			else if(args[i].equals("-labels")){
				labelsPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-outdir")){
				outputPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-config")){
				configPath = args[i+1];
				i++; 
			}
		}
		
		if(countsfilesPath.equals("") || labelsPath.equals("") || outputPath.equals("") || configPath.equals("")){
			System.out.println("Usage Info:\n-countfiles <list of countfiles>\n-labels <label file>\n-outdir <output directory>\n-config <path to config file>");
		}
		else{
	
			Path countsFilePath = Paths.get(countsfilesPath);
			Path labelsFilePath = Paths.get(labelsPath);
			Path configFilePath = Paths.get(configPath);
			
			
			
			File countsFile = countsFilePath.toFile();
			File labelFile = labelsFilePath.toFile();
			File configFile = configFilePath.toFile();

			try {
				DiffExpEval divexpeval = new DiffExpEval(configFile);
				
				divexpeval.getEnrichmentBrowser(countsFile,labelFile,outputPath);
				
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
}
