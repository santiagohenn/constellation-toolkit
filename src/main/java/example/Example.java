package example;

import constellation.tools.ConstellationCoverageComputer;

public class Example {

    public static void main(String[] args) {
        
        String configurationsPath = "./example/example.properties";
        ConstellationCoverageComputer app = new ConstellationCoverageComputer(configurationsPath);
        app.run();

    }
    
}
