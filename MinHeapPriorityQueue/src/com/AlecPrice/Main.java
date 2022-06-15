package com.AlecPrice;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
public class Main {

    public static void main(String[] args) throws FileNotFoundException {
        // write your code here
        int N; // number of vertices
        int M; // number of directed edges
        int S; // starting vertex
        int D; // destination vertex
        //var list = new LinkedList<String>();
        var graph = new AlmostShortestPathGraph<Integer>();

        //String fileName = "src/Test.txt";
       // var file = new File(fileName);
        Scanner input = new Scanner(System.in);

        while (input.hasNextLine()) {
            N = input.nextInt();
            M = input.nextInt();
            if (M == 0 && N == 0) break;

            S = input.nextInt();
            D = input.nextInt();

            for (int i = 0; i < M; i++) {
                int U = input.nextInt();
                int V = input.nextInt();
                int P = input.nextInt();
                var vertexU = graph.produceVertex(U); // graph.new Vertex(U);
                var vertexV = graph.produceVertex(V); // graph.new Vertex(V);
                graph.addEdge(vertexU, vertexV, P);
            }
            var sourceVertex = graph.produceVertex(S);
            var destinationVertex = graph.produceVertex(D);
            var almostShortestPathCost = graph.getAlmostShortestPath(sourceVertex, destinationVertex);
            var cost = almostShortestPathCost == Integer.MAX_VALUE ? -1 : almostShortestPathCost;
            System.out.println(cost);
            graph.clear();
        }
    }


}
class HeapNode<K extends Comparable<? super K>, V extends Comparable<? super V>> {

    private final K key; // The key stored in this HeapNode.
    private V value; // The priority value stored in this HeapNode.

    HeapNode(K newKey, V newValue)  {
        key = newKey;
        value = newValue;
    }

    public K getKey() {
        return key;
    }

    public V getValue() {
        return value;
    }

    public void setValue(V newValue) {
        value = newValue;
    }

    public int compareTo(HeapNode<K, V> node) throws IllegalArgumentException {
        // First compare their priority values, and then their keys.
        if (node.getValue().equals(value)) {
            return getKey().compareTo(node.getKey());
        } else
            return getValue().compareTo(node.getValue());
    }

    public boolean equals(HeapNode<K, V> node) {
        return node.getKey().equals(key) && node.getValue().equals(value);
    }
}// end of HeapNode Class

/**
 * Heapify_Up(index) – moves an element located at the specified index upwards in the heap to correctly reposition an element whose value is less than the value of its parent. This condition may result from removing an element or from changing an element’s value. This method is described on pages 60-61 of the text, and pseudocode is provided on page 61.
 * <p>
 * Heapify_Down(index) – moves an element located at the specified index downwards in the heap to correctly reposition an element whose value is greater than the value of either of its children. This condition may result from removing an element or from changing an element’s value. This method is described on pages 62-63 of the text, and pseudocode is provided on page 63.
 * StartHeap(N) – initializes an empty heap that is set up to store at most N elements. This operation takes O(N) time, as it involves initializing the array that will hold the heap.
 * Insert(item, value) – inserts the item, item, with an ordering value, value, into the heap at the end of the array, then uses Heapify_Up to position the item so as to maintain the heap order. If the heap currently has n elements, this takes O(log n) time.
 * FindMin() – identifies the minimum element in the heap, which is located at index 1, but does not remove it. This takes O(1) time.
 * Delete(index) – deletes the element in the specified heap position by moving the item in the last array position to index, then using Heapify_Down to reposition that item. This is implemented in O(log n) time for heaps that have n elements.
 * <p>
 * ExtractMin() – identifies and deletes the element with the minimum key value, located at index 1, from the heap. This is a combination of the preceding two operations, and so it takes O(log n) time.
 * <p>
 * Delete(item) – deletes the element item form the heap. This can be implemented as a call to Delete(Position[item]), which operates in 𝑂(𝑙𝑜𝑔𝑛)
 * time for heaps that have 𝑛 elements provided Position allows the index of v to be returned in 𝑂(1)
 * time.
 * <p>
 * ChangeKey(item, newValue), which changes the key value of element v to newValue. To implement this operation in 𝑂(𝑙𝑜𝑔𝑛)
 * time, we first need to be able to identify the position of element 𝑣 in the array, which we do by using the Position structure. Once we have identified the position of element 𝑣, we change the key and then apply Heapify-up or Heapify-down as appropriate. You will need to use ChangeKey() method in Part II while implementing Dijkstra Algorithm.
 */
// this is the preferred data structure per the documentation and uses a Binary Heap as a PQ
class MinBinaryHeap<K extends Comparable<? super K>, V extends Comparable<? super V>> {
    private HeapNode<K, V>[] Heap; // heap contains the HeapNodes that make up the binary heap.
    private TreeMap<K, Integer> Position; // Position's values holds the index of the given item K in the heap.
    private int size; // The current number of nodes in the heap.

    //Constructor
    MinBinaryHeap(int newSize) {
        StartHeap(newSize);
        size = 0;
    }

    private void StartHeap(int size) {
        Heap = new HeapNode[size + 1]; // +1 to allow position 0 to be a temporary working element.
        Position = new TreeMap<>();
    }

    public void Insert(K item, V value) throws OutOfMemoryError {

        if (Position.containsKey(item)) {// If the item already exists in the heap, update its value.
            ChangeKey(item, value);
        } else if (isFull()) {
            throw new OutOfMemoryError("Error while executing Insert(K, V) in BinaryMinHeap for item \"" + item + "\", value \"" + value + "\": BinaryMinHeap is full!");
        } else {
            size++;
            Heap[size] = new HeapNode<>(item, value);
            Position.put(item, size);
            Heapify_Up(size);
        }
    }

    public boolean isEmpty() {
        return size == 0;
    }

    public boolean isFull() {
        return size == Heap.length - 1;
    }

    // Function to return the index of the
// parent node of a given node
    static int parent(int index) {
        return (index) / 2;
    }

    // Function to return the index of the
// left child of the given node
    static int leftChild(int index) {
        return ((2 * index));
    }

    // Function to return the index of the
// right child of the given node
    static int rightChild(int index) {
        return ((2 * index) + 1);
    }

    void Heapify_Up(int index) {
        if (index > 1) {
            int parent = Math.floorDiv(index, 2); // index's parent node
            if (Heap[index].compareTo(Heap[parent]) < 0) {
                swap(index, parent, Heap);
                Position.put(Heap[index].getKey(), index);
                Position.put(Heap[parent].getKey(), parent);
                Heapify_Up(parent);
            }
        }
    }

    void Heapify_Down(int index) {
        if (leftChild(index) <= size) { // Skip if 2*index > size.
            int parent = 0; // The eventual child node to swap with the parent node.
            if (leftChild(index) < size) {
                int left = leftChild(index);
                int right = rightChild(index);
                // If the left child node's priority value, then item, is < the right child node's priority value and item.
                if (Heap[left] == null) {
                    parent = left;
                } else if (Heap[right] == null) {
                    parent = right;
                } else if (Heap[left].compareTo(Heap[right]) < 0) // Bring up the smaller of the two.
                    parent = left;
                else {
                    parent = right;
                }
            } else if (leftChild(index) == size) { // If the index is the last index in the heap.
                parent = leftChild(index);
            }
            if (Heap[parent].compareTo(Heap[index]) < 0) { // parent's item's priority value, then item, is less than index's priority value or item.
                HeapNode<K, V> temp = Heap[index];
                Heap[index] = Heap[parent];
                Heap[parent] = temp;
                Position.put(Heap[index].getKey(), index);
                Position.put(Heap[parent].getKey(), parent);
                Heapify_Down(parent);
            }
        }
    }

    // double check to make sure this is accurately swapping the nodes
    void swap(int index, int Parent, HeapNode heap[]) {
        HeapNode<K, V> temp = heap[index];
        heap[index] = heap[Parent];
        heap[Parent] = temp;
        // System.out.println(Arrays.toString(heap)); // print array passed in
    }

    public void ChangeKey(K item, V newValue) throws IllegalArgumentException, IllegalStateException, NullPointerException {
        if (item == null) {
            throw new IllegalArgumentException("Error while executing ChangeKey(K, V) in BinaryMinHeap: The item parameter is null!");
        } else if (newValue == null) {
            throw new IllegalArgumentException("Error while executing ChangeKey(K, V) in BinaryMinHeap for item \"" + item + "\": The newValue parameter is null!");
        } else if (!isEmpty()) {
            Integer pos = Position.get(item); // Fetch the item's position within the heap.
            if (pos == null) {
                throw new NullPointerException("Error while executing ChangeKey(K, V) in BinaryMinHeap for item \"" + item + "\", newValue \"" + newValue + "\": The item does not exist!");
            }
            // could separate this into a decrease and increase key but this is a short handed version to do both for the change key needed in this assignment
            V oldValue = Heap[pos].getValue();
            Heap[pos].setValue(newValue);

            if (newValue.compareTo(oldValue) < 0) {
                Heapify_Up(pos);
            } else if (newValue.compareTo(oldValue) > 0) {
                Heapify_Down(pos);
            }
        } else {
            throw new IllegalStateException("Error while executing ChangeKey(K, V) in BinaryMinHeap for item \"" + item + "\", newValue \"" + newValue + "\": The heap is empty!");
        }
    }

    public K FindMin() {
        if (isEmpty() != false) {
            return null;
        } else {
            return Heap[1].getKey();
        }
    }

    private void Delete(int index) throws IndexOutOfBoundsException, IllegalStateException {
        if (isEmpty()) {
            throw new IllegalStateException("Error while executing Delete(int) in BinaryMinHeap for index \"" + index + "\": The BinaryMinHeap is empty!");
        }
        if (index < size && index > 0) {
            Position.remove(Heap[index].getKey());
            Heap[index] = Heap[size];
            Heap[size] = null;
            size--;
            Position.put(Heap[index].getKey(), index);
            Heapify_Down(index);
        } else if (index == size) {
            Position.remove(Heap[size].getKey());
            Heap[size] = null;
            size--;
        } else {
            throw new IndexOutOfBoundsException("Error while executing Delete(int) in BinaryMinHeap for index \"" + index + "\": index out of bounds (current size: " + size + ")!");
        }
    }

    //same as poll in lib
    public K ExtractMin() {
        if (isEmpty() != false) {
            return null;
        }
        K result = FindMin();
        Delete(1);
        return result;

    }

    public void Delete(K item) throws IllegalArgumentException, IllegalStateException {
        // First check that the item parameter is not null.
        if (item == null) {
            throw new IllegalArgumentException("Error while executing Delete(K) in BinaryMinHeap: The item parameter is null!");
        } else if (isEmpty()) {
            throw new IllegalStateException("Error while executing Delete(K) in BinaryMinHeap for item \"" + item + "\": The BinaryMinHeap is empty!");
        }
        Delete(Position.get(item));
    }
}

class AlmostShortestPathGraph<T extends Comparable<? super T>> {

    private Set<Edge> shortestPathEdges;
    private Map<T, Vertex> vertices;
    private List<Edge> edges;

    //Constructor to initialize empty graphs to store the indexes/paths
    public AlmostShortestPathGraph() {
        this.vertices = new HashMap<T, Vertex>();
        this.edges = new LinkedList<Edge>();
        this.shortestPathEdges = new HashSet<Edge>();
    }

    public void clear() {
        this.shortestPathEdges.clear();
        this.vertices.clear();
        this.edges.clear();
    }


    public Map<T, Vertex> getVertices() {
        return this.vertices;
    }


    public void addVertex(Vertex u) {
        if (!this.vertices.containsKey(u.getLabel())) {
            this.vertices.put(u.getLabel(), u);
        }
    }

    public void removeVertex(Vertex u) {
        this.vertices.remove(u);
    }

    public List<Edge> getEdges() {
        return this.edges;
    }

    public Vertex produceVertex(T label) {
        if (vertices.containsKey(label)) {
            return vertices.get(label);
        }

        // for (var v : this.vertices) {
        //     if (v.label.equals(label)) {
        //         return v;
        //     }
        // }

        return new Vertex(label);
    }

    public void addEdge(Vertex u, Vertex v, Integer p) {
        this.addVertex(u);
        this.addVertex(v);
        var edge = new Edge(u, v, p);
        this.edges.add(edge);
        u.addEdge(edge);
//        // var u = edge.getU();
//        for (var v: this.vertices) {
//            if (v.equals(u)) {
//                u.addEdge(edge);
//                return;
//            }
//        }
        // edge.getU().addEdge(edge);
    }

    public void addEdge(Edge e) {
        this.edges.add(e);
        var u = e.getU();
        for (var entry : this.vertices.entrySet()) {
            var v = entry.getValue();
            if (v.equals(u)) {
                u.addEdge(e);
                return;
            }
        }
        e.getU().addEdge(e);
    }

    public void removeEdge(Edge e) {
        e.getU().removeEdge(e);
        this.edges.remove(e);
    }

    public void removeEdge(Vertex u, Vertex v) {
        for (var edge : this.edges) {
            if (edge.getU().equals(u) && edge.getV().equals(v)) {
                u.removeEdge(edge);
                this.edges.remove(edge);
            }
        }

    }


    public void removeEdges(Vertex u, Vertex v /* , paths */) {
        List<Vertex> nodeQueue = new LinkedList<>();
        nodeQueue.add(v);
        // Vertex currentNode;

        while (0 < nodeQueue.size()) {
            var currentNode = nodeQueue.get(0);
            nodeQueue.remove(0);
            if (!currentNode.equals(u)) {
                for (var ancestor : currentNode.getAncestors()) {
                    var e = ancestor.getEdge(currentNode);
                    ancestor.removeEdge(e);
                    this.edges.remove(e);
                    nodeQueue.add(ancestor);
                }
            }
        }
    }

    public void dijkstra(Vertex source) {
        source.setCost(0);
        var pq = new MinBinaryHeap<Vertex, Integer>(this.vertices.size());
        for (var entry : this.vertices.entrySet()) {
            var v = entry.getValue();
            if (!source.equals(v)) {
                v.setCost(Integer.MAX_VALUE);
                v.clearAncestors();
            }
            pq.Insert(v, v.getCost());
        }

        while (!pq.isEmpty()) {
            var u = pq.ExtractMin();
            for (var v : u.getNeighbors()) {
                var edge = u.getEdge(v);
                var cost = u.getCost();
                var alt = cost == Integer.MAX_VALUE ? cost : cost + edge.getWeight();
                if (alt < v.getCost()) {
                    v.clearAncestors();
                    v.setCost(alt);
                    v.addAncestor(u);
                    pq.ChangeKey(v, alt);
                } else if (alt == v.getCost()) {
                    v.addAncestor(u);
                }
            }
        }
    }

    public Integer getShortestPath(Vertex source, Vertex destination) {
        dijkstra(source);
        var cost = destination.getCost();

        Queue<Vertex> q = new LinkedList<>();
        var v = destination;
        q.offer(v);
        while (0 < q.size()) {
            v = q.poll();
            for (var u : v.getAncestors()) {
                q.offer(u);
                shortestPathEdges.add(u.getEdges().get(v));
            }
        }

        return cost;
    }

    public Integer getAlmostShortestPath(Vertex source, Vertex destination) {
        var optimalCost = getShortestPath(source, destination);

        for (var e : this.shortestPathEdges) {
            removeEdge(e);
        }

        return getShortestPath(source, destination);
    }

    @Override
    public String toString() {
        return "AlmostShortestPathGraph{" +
                "vertices=" + vertices +
                ", edges=" + edges +
                ", shortestPathEdges=" + shortestPathEdges +
                '}';
    }

    class Vertex implements Comparable<Vertex> {
        private T label;
        private Map<Vertex, Edge> edges;
        private List<Vertex> ancestors;
        private int cost;

        Vertex() {
        }

        Vertex(T label) {
            this.label = label;
            edges = new HashMap<>();
            ancestors = new ArrayList<>();
            cost = Integer.MAX_VALUE;
        }

        public Integer getIntegerLabel() {
            return (Integer) this.label;
        }

        public int hashCode() {
            return (Integer) this.label;
        }

        public void clearAncestors() {
            this.ancestors.clear();
        }

        public boolean equals(Vertex o) {
            return this.label.equals(o.label);
        }

        public T getLabel() {
            return this.label;
        }

        public List<Vertex> getAncestors() {
            return this.ancestors;
        }

        public void addAncestor(Vertex v) {
            this.ancestors.add(v);
        }

        public void removeAncestor(Vertex v) {
            this.ancestors.remove(v);
        }

        public Edge getEdge(Vertex v) {
            return this.edges.get(v);
        }

        public Map<Vertex, Edge> getEdges() {
            return this.edges;
        }

        public List<Vertex> getNeighbors() {
            var neighbors = new LinkedList<Vertex>();
            for (var entry : this.edges.entrySet()) {
                neighbors.add(entry.getValue().getV());
            }
            return neighbors;
        }

        public void setEdges(Map<Vertex, Edge> edges) {
            this.edges = edges;
        }

        public void addEdge(Edge edge) {
            this.edges.put(edge.getV(), edge);
        }

        public void removeEdge(Edge edge) {
            this.edges.remove(edge.getV());
        }

        @Override
        public int compareTo(Vertex o) {
            return this.label.compareTo(o.label);
        }

        @Override
        public String toString() {
            return "Vertex: " +
                    label +
                    ' ';
        }

        public void setCost(int cost) {
            this.cost = cost;
        }

        public int getCost() {
            return this.cost;
        }
    }

    class Edge {
        Vertex u;
        Vertex v;
        int weight;

        Edge(Vertex u, Vertex v, int weight) {
            this.u = u;
            this.v = v;
            this.weight = weight;
        }

        @Override
        public String toString() {
            return "Edge: {" +
                    " u = " + u +
                    ", v = " + v +
                    ", weight = " + weight +
                    '}';
        }

        public void setU(Vertex u) {
            this.u = u;
        }

        public Vertex getV() {
            return v;
        }

        public Vertex getU() {
            return u;
        }

        public void setV(Vertex v) {
            this.v = v;
        }

        public int getWeight() {
            return weight;
        }

        public void setWeight(int weight) {
            this.weight = weight;
        }
    }
}