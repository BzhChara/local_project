import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.JFrame;
import javax.swing.JPanel;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.BasicStroke;

class Stock {
    String Name;
    ArrayList<String[]> Data; // Date, Open, High, Low, Close
    Stock(String Name, ArrayList<String[]> Data) {
        this.Name = Name;
        this.Data = Data;
    }

    @Override
    public String toString() {
        StringBuffer str = new StringBuffer();
        for (String[] i:Data) {
            str.append(String.join(" ", i)).append("\n");
        }
        return "Name:" + Name + "\n" + "Data:" + str.toString();
    }
}

class Draw extends JPanel {
    ArrayList<String[]> Data = new ArrayList<>();
    Draw(ArrayList<String[]> Data) {
        this.Data = Data;
        this.setPreferredSize(new Dimension(1280, 1024));
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);

        Graphics2D g2d = (Graphics2D) g;

        setBackground(Color.white);

        BasicStroke thinStroke = new BasicStroke(2);    
        BasicStroke mediumStroke = new BasicStroke(6);   
        
        g2d.setColor(Color.black);
        g2d.drawLine(10, 10, 10, 1014); 
        g2d.drawLine(10, 1014, 1270, 1014); 

        int x = 20;
        for (int i = 0; i < Data.size(); i++) {
            double open = Double.parseDouble( Data.get(i)[1]);
            double high = Double.parseDouble( Data.get(i)[2]);
            double low = Double.parseDouble( Data.get(i)[3]);
            double close = Double.parseDouble( Data.get(i)[4]);

            int openY = (int)(1024 - open * 10) - 10;
            int highY = (int)(1024 - high * 10) - 10;
            int lowY = (int)(1024 - low * 10) - 10;
            int closeY = (int)(1024 - close * 10) - 10;

            if (close > open) {
                g2d.setColor(Color.red);
                g2d.setStroke(mediumStroke);
                
                g2d.drawLine(x, closeY, x, openY);

                g2d.setStroke(thinStroke);
                g2d.drawLine(x, highY, x, lowY);
            }
            else if (close < open) {
                g2d.setColor(Color.green);
                g2d.setStroke(mediumStroke);
                
                g2d.drawLine(x, openY, x, closeY);

                g2d.setStroke(thinStroke);
                g2d.drawLine(x, highY, x, lowY);
            }
            
            else {
                g2d.setColor(Color.black);
                g2d.setStroke(thinStroke);
                g2d.drawLine(x, highY, x, lowY);
                
                g2d.drawLine(x - 2, openY, x + 2, openY);
            }
            g2d.setPaint(Color.black);
            g2d.setFont(new Font("Arial", Font.PLAIN, 12));
            g2d.drawString(Data.get(i)[0], x - 5, 1024);

            x = x + 20;
        }
    }
}

public class Stock_Images {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        ArrayList<Stock> stocks = new ArrayList<>();

        System.out.println("Please enter the number of stocks:");

        int epoch;
        epoch = sc.nextInt();
        sc.nextLine(); // '/n'

        for (int i = 0; i < epoch; i++) {
            System.out.println("Please enter the stock name:");
            
            String name;
            name = sc.nextLine();

            System.out.println("Please enter the relevant data for the corresponding date of this stock(enter 0 to exit):");
            System.out.println("Date Open High Low Close"); // eg:2.30 1 1 1 1 =)

            ArrayList<String[]> data = new ArrayList<>();
            while (true) {
                String tmp;
                tmp = sc.nextLine();
                if (tmp.equals("0")) {
                    break;
                }

                data.add(tmp.split(" "));
            }
            stocks.add(new Stock(name, data));
        }

        System.out.println("----------Stock Information----------");
        
        for (Stock i:stocks) {
            System.out.println(i.toString());

            System.out.print("Generating K-line chart...");

            String name = i.Name;
            ArrayList<String[]> data = i.Data;

            JFrame frame = new JFrame("Stock:" + name);  
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

            Draw drawingStockImage = new Draw(data);

            frame.add(drawingStockImage);
            frame.pack();
            frame.setVisible(true);

            System.out.println("--------------------------------------------------");
        }

        sc.close();
    }
}