class Concurrency:
  def __init__(self, biochip, source, sink, path):
    """
    initializes a flat representation of the chip,
    constructs the ILP model and awaits execution.
    biochip =>  The 3D biochip instance.
    source => A tuple (a,b) of the location of the source
    sink => A tuple (i,j) of the location of the sink
    path => A list of tuples [(a,b), (c,d) ,..., (n,m)]
            that denote off-limit cells.
    """
    x,y = biochip.chip.shape[:-1]
    t = len(path)
    m,n = sink
    xst,yst = source
    #Create biochip as array of dictionaries
    chip = pulp.LpVariable.dicts("ArraySpot",
                      [(i,j,k) for i in range(1,x+1)
                      for j in range(1,y+1)
                      for k in range(1,t+1)],
                      0,1, cat="Integer")
    #Add sink
    for i in range(1,t+1):
      chip[(n,m,i)] = pulp.LpVariable("Sink", 0,1, cat="Integer")

    #Add virtual cell
    for i in range(1, t+1):
      chip[(n+1,m, i)] = pulp.LpVariable("VirtualSpot", 0, 1, cat="Integer")

    #Create problem instance
    self.ilp = pulp.LpProblem('Concurrency Testing ILP', pulp.LpMinimize)

    #OBJECTIVE FUNCTION
    self.ilp += pulp.lpSum(chip[m,n,k]*k for k in range(1,t+1))

    #CONSTRAINTS -- Testing Requirement
    #Any cell (i,j) in the array available for testing
    #should be visited at least once
    self.ilp += pulp.lpSum(chip[i,j,k] for i in range(1,x+1)
                                       for j in range(1,y+1)
                                       for k in range(1,t+1)) >= 1
 
    #Any cell (i,j) in the array that is running biomedical
    #assays cannot be visited by the test droplet
    self.ilp += pulp.lpSum(chip[i,j,k] for i in range(1,x+1)
                                       for j in range(1,y+1)
                                       for k in range(1,t+1)
                                       if (i,j) in path) == 0
 
    #The sink (n,m) should be visited by the test droplet
    #exactly once
    self.ilp += pulp.lpSum(chip[m,n,k] for k in range(1,t+1)) == 1

    #CONSTRAINTS -- Resource Constraint
    #
    self.ilp += pulp.lpSum(chip[i,j,k] for i in range(1, m+1)
                                   for j in range(1,n+1)
                                   for k in range(1,t+1)
                                   if (i,j,k) in chip.keys()) == 1

    #
    for k in range(2, t+1):
      self.ilp += chip[m+1,n,k] == pulp.lpSum(chip[m,n,x] for x in range(1, k-1))

    #CONSTRAINTS -- Starting Point
    #Source should be visited by the test droplet at t == 1
    self.ilp += chip[yst+1,yst,1] == 1

    #CONSTRAINTS == Movement Rules
    #
    rowSum = dict()
    for x in range(1, t+1):
      rowSum[x] = pulp.lpSum(chip[i,j,x]*i for i in range(1,m)
                                      for j in range(1,n))

    colSum = dict()
    for x in range(1, t+1):
      colSum[x] = pulp.lpSum(chip[i,j,k]*j for i in range(1,m)
                                      for j in range(1,n))

    deltaColSum = dict()
    deltaRowSum = dict()
    for x in range(1,t):
      deltaColSum[x] = colSum[x+1] - colSum[x] #FIXME: Absolute Value
      deltaRowSum[x] = rowSum[x+1] - rowSum[x] #FIXME: Absolute Value

    self.ilp += deltaColSum[x] + deltaRowSum[x] <= 1
    return

  def absolute(self, number):
    """
    An attempted workaround for the pulp/abs() clash
    """
    if number < 0:
      return -number
    else:
      return number
  
  def solve(self):
    """
    solves the ILP problem
    """
    return self.ilp.solve()
