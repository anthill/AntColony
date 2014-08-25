"use strict"

var shell = require("game-shell")();
var Map = require("harmony-collections").Map;
var Set = require("harmony-collections").Set;
var dt = require("delaunay-triangulate");

var sqrt = Math.sqrt;
var pow = Math.pow;
var floor = Math.floor;
var random = Math.random;


var nbPoints = 100;
var nbAnts = 200;
var nbCities = 10;

// generate placed points
var cityList = [];
while(cityList.length < nbCities){
      var randomnumber=Math.ceil(Math.random()*nbPoints)
      var found=false;
      for(var i=0;i<cityList.length;i++){
            if(cityList[i]==randomnumber){found=true;break}
      }
      if(!found)cityList[cityList.length]=randomnumber;
}
var citySet = new Set(cityList);


function sign(x) { return x ? x < 0 ? -1 : 1 : 0; }

//Initialize triangulation
var points = new Array(nbPoints)
    for(var i=0; i<nbPoints; ++i) {
        points[i] = {id : i, x:random(), y:random()};
    }
var cells = dt(points.map(function(p){return [p.x, p.y]}))

// create edges
var nextEdges = new Map();
var edges = [];
var permutations = [[0,1], [1,0], [0,2], [2,0], [1,2], [2,1]];
var nbEdges = 0;
cells.forEach(function(cell){
  for(var i=0; i<permutations.length; ++i){
    var s = permutations[i][0];
    var d = permutations[i][1];
    var ps = points[cell[s]];
    var pd = points[cell[d]];
    var edge = {id : nbEdges,
          source: cell[s], 
          destination: cell[d], 
          distance : sqrt( pow(ps.x - pd.x, 2) + pow(ps.y - pd.y, 2) ),
          direction : Math.atan((pd.y-ps.y)/(pd.x-ps.x)),
          orientation : sign((pd.x-ps.x)),
          pheromon : 0
          };
    var nexts;
    if(nextEdges.has(ps.id)){
        nexts = nextEdges.get(ps.id);
        nexts.push(edge);
    } else {
        nexts = [edge];
    }
    nextEdges.set(ps.id, nexts);
    edges.push(edge);
    nbEdges++;
  }
  
})

// initialize ants
var population = new Array(nbAnts);
var i,j;
for (i = 0; i < nbAnts; i++) {
    // take a random edge
    var edge = edges[Math.floor(edges.length*random())];
    var x = points[edge.source].x 
    var y = points[edge.source].y
    population[i] = new Ant(x, y, edge);
}

var canvas, context;

shell.on("init", function() {
    canvas = document.createElement("canvas");
    canvas.width = shell.width;
    canvas.height = shell.height;
    context = canvas.getContext("2d");
    shell.element.appendChild(canvas);
})

shell.on("resize", function(w, h) {
    canvas.width = w;
    canvas.height = h;
})

shell.on("render", function() {
    var w = canvas.width;
    var h = canvas.height;
    var mouse = [shell.mouseX/w, shell.mouseY/h];
    context.setTransform(w, 0, 0, h, 0, 0);
    context.fillStyle = "#fff";
    context.fillRect(0,0,w,h);

    // edges
    context.strokeStyle = "#000";
    for(var i=0; i<edges.length; ++i) {
        var edge = edges[i];
        if (edge.pheromon != 0){
            context.lineWidth = 0.001 * edge.pheromon;
        }else {
            context.lineWidth = 0.00001;
        }
        context.beginPath();
        context.moveTo(points[edge.source].x, points[edge.source].y);
        context.lineTo(points[edge.destination].x, points[edge.destination].y);
        context.stroke();
    }

    // vertices
    for(var i=0; i<points.length; ++i) {
        context.beginPath()
        var point = points[i];
        if (citySet.has(point.id)) {
            context.fillStyle = "#0101DF";
            context.arc(point.x, point.y, 0.006, 0, 2*Math.PI);
        }
        else {
            context.fillStyle = "#000";
            context.arc(points[i].x, points[i].y, 0.003, 0, 2*Math.PI);
        }
        context.closePath();
        context.fill();
    }

    // move ants
    for (i = 0; i < nbAnts; i++) {
        population[i].transit();
    }

    // pheromon evaporation
    for (i = 0; i < edges.length; i++) {
        if(edges[i].pheromon > 0){
            edges[i].pheromon -= 0.001;
        }
    }


    for(var i=0; i<population.length; ++i) {
        context.beginPath()
        var x = population[i].posX //+ 0.01*random();
        var y = population[i].posY //+ 0.01*random();
        if (population[i].state === "pheromoning"){context.fillStyle = "#FF0000"}
        else {context.fillStyle = "#610B0B"}
        context.arc(x, y, 0.003, 0, 2*Math.PI)
        context.closePath()
        context.fill()
    }
  
})

function Ant(x, y, edge) {                                            
    this.posX = x;                
    this.posY = y;
    this.edge = edge;
    this.step = 0;
    this.state = "forage";
    this.transit = statemachine; 
    this.move = move;
    this.edges = [];
    this.lastCity = undefined;
}
// forage: the ant wanders around without any pheromon deposition
// once it finds a city, it starts remembering the nodes it goes through
// when it finds another city, it computes the path length and adds pheromons one each edges
// proportionnaly to the shortestness of the path
// it resets the list of nodes and continues
// while foraging the ant choses the path with a pheromon preference

function move() {
    var edgeChanged;
    var cityReached = false;
    // on edge
    if (this.step < this.edge.distance){
        this.posX += 0.005*Math.cos(this.edge.direction)*this.edge.orientation;
        this.posY += 0.005*Math.sin(this.edge.direction)*this.edge.orientation;
        this.step += 0.005;
        edgeChanged = false;
    // on vertex
    } else {
        this.step = 0;
        this.posX = points[this.edge.destination].x;
        this.posY = points[this.edge.destination].y;
        var possibleEdges = nextEdges.get(this.edge.destination);
        // flip a coin and either take the smelliest path of a random one
        if (random() > 0.5){
            var smells = possibleEdges.map(function(e){return e.pheromon});
            var index = smells.indexOf(Math.max.apply(Math, smells));
            this.edge = possibleEdges[index];
        } else {
            this.edge = possibleEdges[floor(random()*possibleEdges.length)];
        }
        cityReached = citySet.has(this.edge.source);
        edgeChanged = true;
    }
    return {cityReached: cityReached, edgeChanged: edgeChanged};
}


function statemachine() {
    switch (this.state) {
        case "forage":
            var res = this.move();
            if (res.cityReached) {
                this.state = "pheromoning";
                this.lastCity = this.edge.source;
            };
            break;
        case "pheromoning":
            var res = this.move();
            if (res.edgeChanged) {
                this.edges.push(this.edge);
                // found a city
                if (res.cityReached && (this.edge.source != this.lastCity) ){
                    // compute the length of the path
                    var pathLength = this.edges.map(function(e){return e.distance}).reduce(function(a,b){return a + b});
                    var deltaPheromone = 1/pathLength;
                    this.edges.forEach(function(e){e.pheromon += deltaPheromone});
                    // console.log(deltaPheromone, this.edges);
                    this.edges = [this.edge];
                    this.lastCity = this.edge.source;
                }
            }
          break;

    }
}