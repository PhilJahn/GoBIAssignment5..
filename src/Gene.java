
public class Gene {

	private String name;
	
	private double fc;
	
	public Gene(String name){
		this.name = name;
	}
	
	
	public String getName(){
		return name;
	}
	
	public void setFC (double fc){
		this.fc = fc;
	}
	
	public double getFC(){
		return fc;
	}
}
