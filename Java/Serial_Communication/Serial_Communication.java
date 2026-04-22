package Serial_Communication;

public class Serial_Communication {
    private boolean running = true;
    private final StringBuilder buffer = new StringBuilder();

    public void start() {
        try {
            openSerialPort();

            while (running) {
                String data = readFromSerial();

                if (data != null && !data.isEmpty()) {
                    buffer.append(data);

                    String packet;
                    while ((packet = extractPacket()) != null) {
                        Integer result = parsePacket(packet);
                        System.out.println("解析结果: " + result);
                    }
                }
            }
        } catch (Exception e) {
            System.out.println("发生异常: " + e.getMessage());
            e.printStackTrace();
        } finally {
            closeSerialPort();
        }
    }

    private void openSerialPort() {
        System.out.println("打开串口");
    }

    private String readFromSerial() {
        // 串口读取逻辑
        return null;
    }

    private String extractPacket() {
        // 按需修改
        int index = buffer.indexOf("\n");
        if (index == -1) {
            return null;
        }

        String packet = buffer.substring(0, index).trim();
        buffer.delete(0, index + 1);
        return packet;
    }

    private Integer parsePacket(String packet) {
        // 按需修改
        return packet.length();
    }

    private void closeSerialPort() {
        System.out.println("关闭串口");
    }

    public static void main(String[] args) {
        new SerialCommunication().start();
    }
}
