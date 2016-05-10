import java.util.*;
import java.io.*;


class random{
    public static void main(String[] args) throws FileNotFoundException, IOException{
        int size = args.length == 0 ? 10 : Integer.parseInt(args[0]);
        boolean[][] matrix = new boolean[size][size];

        int c = (int)(Math.random() * 5);
        for(int i = 0; i < size * c; i++){
            int u = (int)(Math.random() * size);
            int v = (int)(Math.random() * size);
            if(u == v) continue;
            matrix[u][v] = true;
            matrix[v][u] = true;
        }

        
        PrintWriter w = new PrintWriter(new BufferedWriter(new FileWriter("graph.txt")),true);
        w.println(size);
        for(int i = 0; i < size; i++){
            for(int j = 0; j < size; j++){
                w.print(matrix[i][j] ? "1 " : "0 ");
            }
            w.println();
        }

    }
}
