import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// The problem of duplicate input was not handled.

class Data implements Serializable {
    String Name;
    long ID;
    String Designation;
    int Experience;
    int Age;
    Data(String Name, long ID, String Designation, int Experience, int Age) {
        this.Name = Name;
        this.ID = ID;
        this.Designation = Designation;
        this.Experience = Experience;
        this.Age = Age;
    }

    @Override  
    public String toString () {
        return "Name:" + Name + "\n" + "ID:" + ID + "\n" + "Designation:" + Designation + "\n" +  "Experience:" + Experience + "\n" + "Age:" + Age;
    }
}   

class Function {
    // All employees
    public void Show (ArrayList<Data> employee) {
        for (int i = 0;i < employee.size(); i++) {
            System.out.println("\tEmployee" + (i+1) + ":\t");
            System.out.println(employee.get(i).toString());
        }
    }
    // Single employee
    public void Show (Data employee) {
        System.out.println(employee.toString());
    }

    public void Insert (ArrayList<Data> employee, Scanner sc) {
        String st;

        System.out.print("Input the new data:");
        st = sc.nextLine();

        Pattern p = Pattern.compile("^\\s*\\w+\\s+\\d+\\s+\\w+\\s+\\d+\\s+\\d+\\s*$");
        Matcher m = p.matcher(st);
        while (!m.matches()) {
            System.out.print("The input format is incorrect, please re-enter:");
            st = sc.nextLine();
            m = p.matcher(st);
        }

        employee.add(new Data(st.split("\\s+")[0], Long.parseLong(st.split("\\s+")[1]), st.split("\\s+")[2], Integer.parseInt(st.split("\\s+")[3]), Integer.parseInt(st.split("\\s+")[4])));
    }

    public void Delete (ArrayList<Data> employee, Scanner sc) {
        String st;

        System.out.print("Input the name of the person whose information you want to delete:");
        st = sc.nextLine();

        Pattern p = Pattern.compile("\\s*\\w+\\s*");
        Matcher m = p.matcher(st);

        while (!m.matches()) {
            System.out.print("The input format is incorrect, please re-enter:");
            st = sc.nextLine();
            m = p.matcher(st);
        }

        boolean bool = true;
        for (Data tmp:employee) {
            if (tmp.Name.equals(st)) {
                employee.remove(tmp);
                System.out.println("Delete successfully!");
                bool = false;
                break;
            }
        }

        if (bool) {
            System.out.println("Target not found, please check your input!");
        }
    }

    public void Search (ArrayList<Data> employee, Scanner sc) {
        String st;

        System.out.print("Input the name of the person whose information you want to search:");
        st = sc.nextLine();

        Pattern p = Pattern.compile("\\s*\\w+\\s*");
        Matcher m = p.matcher(st);

        while (!m.matches()) {
            System.out.print("The input format is incorrect, please re-enter:");
            st = sc.nextLine();
            m = p.matcher(st);
        }

        boolean bool = true;
        for (Data tmp:employee) {
            if (tmp.Name.equals(st)) {
                Show(tmp);
                bool = false;
                break;
            }
        }

        if (bool) {
            System.out.println("Target not found, please check your input!");
        }
    }
}

public class Program_for_Employee_Management_System {
    @SuppressWarnings("unchecked")
    public static void main(String[] args) {
        ArrayList<Data> employee = new ArrayList<>();

        Scanner sc = new Scanner(System.in);

        // Load data
        File file = new File("employee_data.ser");
        if (file.exists()) {
            try (FileInputStream fis = new FileInputStream(file); ObjectInputStream ois = new ObjectInputStream(fis)) {
                Object loadedData = ois.readObject();
                if (loadedData instanceof ArrayList) {
                    employee = (ArrayList<Data>) loadedData;
                    System.out.println("Loaded " + employee.size() + " employees from file.");
                }
                
            } catch (FileNotFoundException e) {
                System.err.println("File not found: " + e.getMessage());
            } catch (IOException e) {
                System.err.println("Error reading file: " + e.getMessage());
            } catch (ClassNotFoundException e) {
                System.err.println("Class not found during deserialization: " + e.getMessage());
            }
        }

        if (employee.isEmpty()) {
            System.out.println("----------Built Employee Table----------");
            System.out.print("Input the number of employees:");

            int epoch;
            String st;
            epoch = sc.nextInt();
            st = sc.nextLine();

            System.out.println("Input the property of employees(Name, ID, Designation, Experience, Age):");
            
            for (int i = 0; i < epoch; i++) {
                System.out.print("Employee" + (i+1) + ":");
                st = sc.nextLine();

                Pattern p = Pattern.compile("^\\s*\\w+\\s+\\d+\\s+\\w+\\s+\\d+\\s+\\d+\\s*$");
                Matcher m = p.matcher(st);
                while (!m.matches()) {
                    System.out.print("The input format is incorrect, please re-enter:");
                    st = sc.nextLine();
                    m = p.matcher(st);
                }

                employee.add(new Data(st.split("\\s+")[0], Long.parseLong(st.split("\\s+")[1]), st.split("\\s+")[2], Integer.parseInt(st.split("\\s+")[3]), Integer.parseInt(st.split("\\s+")[4])));
            }
        }

        System.out.println("----------Menu----------");
        System.out.println("0:Exit\n1:Show Employee Table\n2:Insert New Entries\n3:Delete The Entries\n4:Search A Record");

        int num = -1;
        Function fun = new Function();
        
        while (num != 0) {
            System.out.print("Input number selection function:");
            num = sc.nextInt();
            sc.nextLine();

            if (num == 1) {
                System.out.println("----------Employee Table----------");
                fun.Show(employee);
            }
            else if (num == 2) {
                fun.Insert(employee, sc);
            }
            else if (num == 3) {
                fun.Delete(employee, sc);
            }
            else if (num == 4) {
                fun.Search(employee, sc);
            }
        }

        // Save data
        try (FileOutputStream fos = new FileOutputStream(file); ObjectOutputStream oos = new ObjectOutputStream(fos)) {
            oos.writeObject(employee);
            System.out.println("Employee data saved successfully!");  
        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + e.getMessage());
        } catch (IOException e) {
            System.err.println("Error saving data: " + e.getMessage());
        }

        System.out.println("Successfully exited the system!");
        sc.close();
    }
}