package constellation.tools.utilities;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TimeZone;

public class TimeUtils {

    public static long stamp2unix(String timestamp) {
        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSS");
        dateFormat.setTimeZone(TimeZone.getTimeZone("UTCG"));
        Date parsedDate = new Date();

        try {
            parsedDate = dateFormat.parse(timestamp);
            return parsedDate.getTime();
        } catch (ParseException parseException) {
            parseException.printStackTrace();
            return parsedDate.getTime();
        }
    }

}
