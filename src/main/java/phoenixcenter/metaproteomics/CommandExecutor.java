package phoenixcenter.metaproteomics;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.Setter;
import lombok.extern.log4j.Log4j2;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;

@Log4j2
@AllArgsConstructor
public class CommandExecutor {

    private ExecutorService threadPool;

    /**
     * Execute system commands.
     * <p>
     * Note: do not merge error and output, we need obtain output not error to do some process sometimes.
     *
     * @param command
     * @param outConsumer the output of command
     * @param errConsumer the error output of command
     */
    public void exec(String command,
                     Consumer<String> outConsumer,
                     Consumer<String> errConsumer) {
        log.info("command>>{}", command);
        try {
            final Process process = Runtime.getRuntime().exec(command);
            Runnable stdOut = () -> {
                try (BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                    String line;
                    if (outConsumer != null) {
                        while ((line = br.readLine()) != null) {
                            outConsumer.accept(line);
                        }
                    } else {
                        while ((line = br.readLine()) != null) {
                            System.out.println(line);
                        }
                    }
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            };

            Runnable stdErr = () -> {
                try (BufferedReader br = new BufferedReader(new InputStreamReader(process.getErrorStream()))) {
                    String line;
                    if (errConsumer != null) {
                        while ((line = br.readLine()) != null) {
                            errConsumer.accept(line);
                        }
                    } else {
                        while ((line = br.readLine()) != null) {
                            System.err.println(line);
                        }
                    }
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            };
            threadPool.submit(stdOut);
            threadPool.submit(stdErr);
            if (process.waitFor() != 0) {
                throw new RuntimeException("command <" + command + "> cannot work");
            }
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * @param command
     * @param outConsumer
     * @see #exec(String, Consumer, Consumer)
     */
    public void exec(String command, Consumer<String> outConsumer) {
        exec(command, outConsumer, null);
    }

    /**
     * @param command
     * @see #exec(String, Consumer, Consumer)
     */
    public void exec(String command) {
        exec(command, null, null);
    }
}
