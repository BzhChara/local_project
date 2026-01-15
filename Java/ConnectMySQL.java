import java.sql.*;

public class ConnectMySQL {
    public static void main(String[] args) {
        String url = "jdbc:mysql://localhost:3306/runoob";
        String username = "root";
        String password = "112303@qq.comV2";
        
        try {
            // 1. 加载驱动
            Class.forName("com.mysql.cj.jdbc.Driver");
            System.out.println("动加载成功");
            
            // 2. 建立连接
            Connection conn = DriverManager.getConnection(url, username, password);
            System.out.println("连接MySQL成功");
            
            // 3. 测试查询
            Statement stmt = conn.createStatement();
            ResultSet rs = stmt.executeQuery("SELECT 'Hello MySQL' as message");
            
            if (rs.next()) {
                System.out.println("测试查询结果: " + rs.getString("message"));
            }
            
            // 4. 关闭连接
            conn.close();
            System.out.println("连接已关闭");
            
        } catch (ClassNotFoundException e) {
            System.out.println("找不到MySQL驱动");
            System.out.println("请确保mysql-connector-java.jar已添加到classpath");
        } catch (SQLException e) {
            System.out.println("数据库连接失败: " + e.getMessage());
        }
    }
}