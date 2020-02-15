package phoenixcenter.metaproteomics;


import lombok.extern.log4j.Log4j2;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Properties;

@Log4j2
public class GlobalConfig {

    private static Properties probs;

    static {
        try {
            init();
        } catch (IOException e) {
            throw new IllegalStateException("Global setting initialization error!!");
        }
    }

    /**
     * Initialize Properties
     *
     * @throws IOException
     */
    public static void init() throws IOException {
        probs = new Properties();
        probs.load(GlobalConfig.class.getResourceAsStream("/mpq.properties"));
    }

    public static String getValue(String key) {
        return probs.getProperty(key);
    }

    public static int getIntValue(String key) {
        return Integer.parseInt(probs.getProperty(key));
    }
}
