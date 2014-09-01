"use strict"

var shell = require("game-shell")();
var Map = require("harmony-collections").Map;
var Set = require("harmony-collections").Set;
var dt = require("delaunay-triangulate");
var parse = require('parse-svg-path')


var nbRandomPoints = 300;
var nbAnts = 1000;
var antVelocity = 0.005;
var textMesh = true;
var nbCity = 10;
var nbStartPoints = 50;

var sqrt = Math.sqrt;
var pow = Math.pow;
var floor = Math.floor;
var random = Math.random;

// Frame definition
var xInit = 0.15, yInit = 0.15;
var w = 0.6,
    h = 0.6;

var svgString = "m 1028.3224,362.79152 c -16.6587,-2.34387 -33.35995,-6.23515 -49.25137,-12.08831 -15.89142,-5.85316 -30.97305,-13.6682 -44.3926,-23.85957 7.4743,-16.05278 16.894,-31.06927 24.4232,-47.05212 17.65805,13.38324 38.61329,24.26005 60.57897,29.97765 21.9657,5.7176 44.9419,6.27598 66.6418,-0.97765 14.2804,-6.92169 17.3033,-19.64436 13.4946,-31.15028 -3.8087,-11.50591 -14.4489,-21.79507 -27.4946,-23.84972 -10.162,-4.07686 -21.1052,-7.14863 -32.1975,-10.11097 -11.0924,-2.96235 -22.334,-5.81526 -33.0931,-9.45439 -10.7591,-3.63914 -21.03562,-8.06449 -30.19779,-14.17171 -9.16216,-6.10723 -17.20996,-13.89632 -23.51158,-24.26293 -5.1794,-10.97058 -7.48544,-22.77583 -7.30512,-34.52966 0.18031,-11.75382 2.84698,-23.4562 7.61299,-34.22103 4.76601,-10.76483 11.63137,-20.5921 20.20905,-28.595692 8.57769,-8.003592 18.86771,-14.183506 30.48305,-17.653621 12.1821,-4.024145 24.8777,-6.357009 37.6988,-7.121901 12.8211,-0.764892 25.7677,0.03819 38.4517,2.285933 12.6841,2.247744 25.1055,5.940153 36.8765,10.953918 11.7709,5.013764 22.8912,11.348885 32.973,18.882053 -8.949,14.00735 -15.1345,30.81362 -25.2587,43.41697 -15.0728,-11.15158 -33.5009,-20.54491 -52.6114,-24.91726 -19.1105,-4.37235 -38.9034,-3.72372 -56.7058,5.20861 -10.0321,7.02789 -13.3348,18.84695 -11.1561,29.53596 2.1787,10.68902 9.8387,20.24799 21.732,22.75572 10.1335,4.43307 21.0241,7.69407 32.1101,10.69023 11.0861,2.99616 22.3676,5.72747 33.283,9.10117 10.9155,3.3737 21.4648,7.38978 31.0865,12.95547 9.6217,5.56569 18.3157,12.68099 25.5204,22.25313 6.6301,9.71635 10.5587,20.88023 12.0064,32.41231 1.4476,11.53207 0.4142,23.43234 -2.8797,34.62146 -3.2939,11.18912 -8.8484,21.66709 -16.443,30.35457 -7.5946,8.68749 -17.2293,15.58449 -28.6837,19.61166 -26.0183,11.85857 -56.0424,10.83764 -84,9 z m -927.99998,-4 c 9.2022,-23.42094 18.429,-46.83835 27.66563,-70.25398 9.23663,-23.41563 18.48308,-46.82948 27.72457,-70.2433 9.24149,-23.41382 18.47802,-46.82761 27.69481,-70.2431 9.21679,-23.41549 18.41385,-46.83269 27.57639,-70.253343 25.99846,-6.71571 45.20583,-13.853157 57.3386,25.993723 8.35197,21.24228 16.73731,42.47074 25.12749,63.69719 8.39019,21.22645 16.78522,42.45089 25.15658,63.68514 8.37136,21.23425 16.71905,42.47831 25.01455,63.74398 8.29549,21.26567 16.53879,42.55296 24.70138,63.87369 -5.90896,6.95561 -24.61617,1.11298 -35.59372,3 -16.75248,2.84859 -26.96421,-0.41416 -28.40628,-19 -6.35452,-16.84315 -12.71521,-33.83068 -18.63059,-50.7668 -16.46723,-0.80053 -35.17768,-2.74281 -53.33727,-3.11439 -18.1596,-0.37158 -35.76834,0.82753 -50.03214,6.30976 -4.15679,10.93124 -8.13806,21.9269 -12.15785,32.90546 -4.0198,10.97855 -8.07813,21.94001 -12.38903,32.80282 -19.05516,1.65495 -38.30798,0.54702 -57.45312,0.86315 l 0,-2.19791 z m 174,-109 c -3.17462,-9.51136 -6.44835,-18.99174 -9.7543,-28.46224 -6.81559,-19.51412 -6.24202,-25.37946 -12.3174,-43.60788 -3.16722,-9.51543 -13.60257,-32.32866 -16.49973,-41.92988 -6.08989,12.88576 -11.132,26.59272 -15.93415,40.47476 -4.80215,13.88203 -9.36435,27.93915 -14.49442,41.52524 -2.39425,12.05801 -22.8815,40.11116 2.18745,34 22.10513,-0.88971 45.90732,2.22185 66.81255,-2 z m 126,-31 c 0,-23.83334 0,-47.66667 0,-71.5 0,-23.83333 0,-47.666667 0,-71.500003 14.64729,1.05313 31.03193,-2.08808 44.60679,1.55415 12.36152,15.221993 24.69446,30.481383 36.98987,45.777973 12.2954,15.29658 24.55326,30.63036 36.76464,46.00112 12.21137,15.37076 24.37625,30.77851 36.48569,46.22304 12.10943,15.44453 24.16342,30.92584 36.15301,46.44372 0.4558,-15.4914 0.72655,-30.98725 0.88119,-46.48582 0.15464,-15.49858 0.19319,-30.99988 0.18458,-46.50218 -0.009,-15.5023 -0.0643,-31.0056 -0.0983,-46.50817 -0.034,-15.50258 -0.0461,-31.004432 0.0325,-46.503833 18.66667,0 37.33333,0 56,0 0,23.833333 0,47.666663 0,71.500003 0,23.83333 0,47.66667 0,71.5 0,23.83333 0,47.66667 0,71.5 0,23.83333 0,47.66667 0,71.5 -14.6441,-1.05348 -31.026,2.08847 -44.59764,-1.55416 -12.54493,-14.76073 -24.84759,-29.74621 -37.03234,-44.83972 -12.18476,-15.0935 -24.2516,-30.29502 -36.32488,-45.48783 -12.07329,-15.1928 -24.15301,-30.37688 -36.36353,-45.43551 -12.21052,-15.05863 -24.55183,-29.9918 -37.14828,-44.68278 -0.22791,15.16496 -0.36674,30.3308 -0.4489,45.49718 -0.0822,15.16637 -0.10766,30.33329 -0.10889,45.5004 -0.001,15.16712 0.0218,30.33443 0.0367,45.50162 0.0149,15.16718 0.0216,30.33422 -0.0122,45.5008 -18.66667,0 -37.33333,0 -56,0 0,-23.83334 0,-47.66667 0,-71.5 0,-23.83333 0,-47.66667 0,-71.5 z m 376,24 c 0,-19.83334 0,-39.66667 0,-59.5 0,-19.83333 0,-39.66666 0,-59.5 -15.33334,0 -30.66667,0 -46,0 -15.33333,0 -30.66667,0 -46,0 0,-16 0,-32.000003 0,-48.000003 19.83333,0 39.66667,0 59.5,0 19.83333,0 39.66667,0 59.5,0 19.83333,0 39.66667,0 59.5,0 19.83333,0 39.66667,0 59.5,0 0,16 0,32.000003 0,48.000003 -15,0 -30,0 -45,0 -15,0 -30,0 -45,0 -0.14707,19.71176 -0.13037,39.43185 -0.0675,59.1545 0.0629,19.72265 0.17199,39.44785 0.20971,59.16984 0.0377,19.72199 0.004,39.44075 -0.21855,59.15053 -0.22261,19.70977 -0.63417,39.41055 -1.35226,59.09656 -14.79248,-1.99928 -46.24135,7.86759 -54.57143,-2.96725 0,-19.1007 0,-38.20139 0,-57.30209 0,-19.1007 0,-38.2014 0,-57.30209 z";

function svgToPoints(svgString) {
    var points = [];
    var edges = [];

    var X = 0;
    var Y = 0;
    var nbPoints = 0;
    var prevPoint;
    var commands = parse(svgString)
    for (var i=0; i<commands.length; i++){
        var command = commands[i];
        switch (command[0]) {
            case "m":
                X += command[1];
                Y += command[2];
                prevPoint = undefined;
                break;
            case "M":
                X = command[1];
                Y = command[2];
                prevPoint = undefined;
                break;  
            case "c":
                X += command[5];
                Y += command[6];
                points.push({id:nbPoints, x:X, y:Y});
                nbPoints++;
                if (prevPoint) {
                    edges.push([prevPoint, nbPoints]);
                    prevPoint = nbPoints;
                }
                break;    
        }
    }
    return {points : points, edges : edges};
}

function range(start, count) {
    return Array.apply(0, Array(count)).map(function (element, index) { return index + start });
}

function sign(x) { return x ? x < 0 ? -1 : 1 : 0; }

// initialize points
var points = [];
var citySet;

if (textMesh){

    var myText = svgToPoints(svgString);
    points = myText.points;
    citySet = new Set(range(0, points.length));
    var scaleX = 0.4;
    var scaleY = 0.2;

    // scale points to [0,1] + scale
    var maxX = Math.max.apply(Math, points.map(function(p){return p.x}));
    var minX = Math.min.apply(Math, points.map(function(p){return p.x}));
    var maxY = Math.max.apply(Math, points.map(function(p){return p.y}));
    var minY = Math.min.apply(Math, points.map(function(p){return p.y}));
    points = points.map(function(p){
        return {
            id:p.id, 
            x: scaleX*(p.x-minX)/(maxX-minX)+0.25, 
            y: scaleY*(p.y-minY)/(maxY-minY)+0.25}
    });

    //add random points
    var nbPoints = points.length;
    for(var i=0; i<nbRandomPoints; ++i) {
        points.push({
            id : nbPoints,
            x:random(),
            y:random()
        });
        nbPoints++;
    }

} else {
    //only add random points
    var nbPoints = 0;
    for(var i=0; i<nbRandomPoints; ++i) {
        points.push({
            id : nbPoints, 
            x:random(), 
            y:random()
        });
        nbPoints++;
    }
    citySet = new Set(range(0, nbCity));
}

// triangulate
var cells = dt(points.map(function(p){
    return [p.x, p.y]
}))

// create edges
var nextEdges = new Map(); // map: Point -> List[Edge] 
var edges = []; // list of all edges
var nbEdges = 0;
cells.forEach(function(cell){
   
    for (var i = 0; i < 3; ++i){  // for each point.id listed in current cell
        var pt = cell[i];
        var nexts = []; // edges accessible by pt

        for (var j = 1; j < 3; ++j){ 

            var ptj = cell[( i + j ) % 3];
            var temp = []; // equivalent of nexts, but for point ptj
            var newEdge = undefined;

            // if key already exists, add the corresponding edges
            if ( nextEdges.has(points[pt]) ){
                var t = nextEdges.get(points[pt]);
                // flatMap to have all possible connected points
                var temppoints = t.map(function(e){ return [e.pt1, e.pt2]})
                                .reduce(function(a, b){return a.concat(b)});

                if (temppoints.indexOf(points[ptj]) == -1){
                    newEdge = createEdge(points[pt], points[ptj]);
                    nexts.push(newEdge);
                    edges.push(newEdge);
                }
            }
            else {
                newEdge = createEdge(points[pt], points[ptj]);
                nexts.push(newEdge);
                edges.push(newEdge);
            }
            if (newEdge != undefined){
                temp.push(newEdge);
                addToNextEdges(points[ptj], cell, temp, nextEdges); // add also the edge to the edge's other point's nextEdges
            }
        }

        // force the edges of the text to be edges of the graph
        if (textMesh && citySet.has(pt)) {
            //console.log('verif');
            var textEdge = createEdge(points[pt], points[pt+1]);
            //var temp2 = [];
            edges.push(textEdge);
            nexts.push(textEdge);
            //temp2.push(textEdge);
            //addToNextEdges(points[pt+1], undefined, temp2, nextEdges);
        }

        //console.log("longueur: " + nexts.length);
        addToNextEdges(points[pt], cell, nexts, nextEdges);
    }
    //console.log(edges.length);
})

//console.log("taille: " + nextEdges.size());

function createEdge(a, b){
    var distance = sqrt( pow(a.x - b.x, 2) + pow(a.y - b.y, 2) );
    var edge = {
        id: nbEdges,
        pt1: a, 
        pt2: b, 
        distance: distance,
        direction: Math.atan((b.y-a.y)/(b.x-a.x)),
        pheromon: 1/distance
        //other: 
    };
    nbEdges ++;
    return edge;
}

function addToNextEdges(a, cell, nexts, nextEdgesMap) {

    if(nextEdgesMap.has(a))
    {
        var existingNextEdges = nextEdgesMap.get(a);
        
        // for every already existing next edge
        for (var i = 0; i < existingNextEdges.length; i++)
        {
            // check if current existing edge is part of current cell
            if (cell != undefined && cell.indexOf(getOtherPoint(existingNextEdges[i], a)) == -1)
                // if not, i.e if current edge hasn't been already added, push current edge
                nexts.push(existingNextEdges[i]);
        }
        
    }

    // update nextEdges with new identified edges
    nextEdgesMap.set(a, nexts);

    //console.log("NextEdges de: " + a.id);
    //console.log(nextEdges.get(a));
}

function getOtherPoint(edge, point){
    // function that returns the other point of an edge
    if (point == edge.pt1)
        return edge.pt2;
    else
        return edge.pt1;
}

// initialize ants
var population = new Array(nbAnts);
var i,j;
var possibleStartPointsId = [];

function initializeStartPoints(){
    for (var i = 0; i < nbStartPoints; i++)
    {
        possibleStartPointsId.push(Math.floor(nbRandomPoints * random()));
    }
}

function generateRandStartPoint(){
    var randId = Math.floor(possibleStartPointsId.length * random());
    var randStartPoint = points[possibleStartPointsId[randId]];
    return randStartPoint;
}

initializeStartPoints();

for (i = 0; i < nbAnts; i++) {
    /*// take a random edge
    var edge = edges[Math.floor(edges.length*random())];
    var x = points[edge.source].x 
    var y = points[edge.source].y*/
    var newAnt = new Ant(generateRandStartPoint());
    newAnt.setDirection();
    population[i] = newAnt;
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
            context.lineWidth = Math.min(0.00001 * edge.pheromon, 0.01);
        } else {
            context.lineWidth = 0.00001;
        }
        context.beginPath();
        context.moveTo(points[edge.pt1.id].x, points[edge.pt1.id].y);
        context.lineTo(points[edge.pt2.id].x, points[edge.pt2.id].y);
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
            edges[i].pheromon -= 0.0001;
        }
    }


    for(var i=0; i<population.length; ++i) {
        context.beginPath()
        var x = population[i].posX + 0.005*random();
        var y = population[i].posY + 0.005*random();
        // var x = population[i].posX;
        // var y = population[i].posY;
        if (population[i].state === "pheromoning"){
            context.fillStyle = "#FF0000";
        }
        else {context.fillStyle = "#610B0B"}
        context.arc(x, y, 0.003, 0, 2*Math.PI)
        context.closePath()
        context.fill()
    }
  
})

function Ant(point) {                                            
    this.posX = point.x;                
    this.posY = point.y;
    this.velocity = antVelocity;
    this.edge = undefined;
    this.step = 0;
    this.state = "forage";
    this.edges = [];
    this.lastCity = undefined;
    this.origin = point;
    this.destination = undefined;
    this.orientation = undefined;
    // methods
    this.transit = statemachine; 
    this.move = move;
    this.setDirection = setDirection;
}
// forage: the ant wanders around without any pheromon deposition
// once it finds a city, it starts remembering the nodes it goes through
// when it finds another city, it computes the path length and adds pheromons one each edges
// proportionnaly to the shortestness of the path
// it resets the list of nodes and continues
// while foraging the ant choses the path with a pheromon preference

function setDirection() {

    var possibleEdges = [];

    for (var i = 0; i < nextEdges.get(this.origin).length; i++)
    {
        possibleEdges[i] = nextEdges.get(this.origin)[i];
    } 

    possibleEdges.splice(possibleEdges.indexOf(this.edge),1);

    // flip a coin and either take the smelliest path of a random one
    if (random() > 0.5){
        var smells = possibleEdges.map(function(e){return e.pheromon});
        var index = smells.indexOf(Math.max.apply(Math, smells));
        this.edge = possibleEdges[index];
    } else {
        this.edge = possibleEdges[floor(random()*possibleEdges.length)];
    }
    // set the destination point, being edge.pt1 or edge.pt2
    this.destination = (this.origin == this.edge.pt1) ? this.edge.pt2 : this.edge.pt1;
    this.orientation = sign((this.destination.x-this.origin.x));
}


function move() {
    var edgeChanged;
    var cityReached = false;
    // on edge
    if (this.step < this.edge.distance){
        this.posX += this.velocity*Math.cos(this.edge.direction)*this.orientation;
        this.posY += this.velocity*Math.sin(this.edge.direction)*this.orientation;
        this.step += this.velocity;
        edgeChanged = false;
    // on vertex
    } else {
        this.step = 0;
        this.origin = this.destination;
        this.posX = this.origin.x;
        this.posY = this.origin.y;

        this.setDirection();

        cityReached = citySet.has(this.origin.id);
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
                this.lastCity = this.origin.id;
            };
            break;
        case "pheromoning":
            var res = this.move();
            if (res.edgeChanged) {
                this.edges.push(this.edge);
                // found a city
                if (res.cityReached && (this.origin.id != this.lastCity) ){
                    // compute the length of the path
                    var pathLength = this.edges.map(function(e){return e.distance}).reduce(function(a,b){return a + b});
                    var deltaPheromone = 1/pathLength;
                    this.edges.forEach(function(e){
                        var a = e.pt1, b = e.pt2, poids = 1;  
                        // increased deposited pheromons for textEdges
                        if (citySet.has(a.id) && citySet.has(b.id) && Math.abs(a.id - b.id) == 1) {
                            poids *= 10;
                        }
                        e.pheromon += (deltaPheromone * poids);
                    });
                    // console.log(deltaPheromone, this.edges);
                    this.edges = [this.edge];
                    this.lastCity = this.origin.id;
                }
            }
          break;

    }
}