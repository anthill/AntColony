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

var svgString = "m 1246.3864,422.67046 c -32.6566,-3.05977 -80.9202,-1.60892 -64,-49.23572 -1.4759,-17.26999 -2.9389,-34.52697 19.6179,-27.76428 28.2051,-7.71782 5.6299,36.53428 29.3821,44 15.6401,11.08789 36.8268,15.16788 55.7601,11.43591 18.9333,-3.73198 35.6132,-15.27593 42.2399,-35.43591 11.114,-30.53447 -5.455,-63.66997 -34.7792,-76 -21.2379,-12.18874 -45.1963,-24.31647 -65.5186,-39.83315 -20.3223,-15.51668 -37.0085,-34.42231 -43.7022,-60.16685 -5.2777,-23.82651 1.2957,-47.48151 15.0269,-66.10238 13.7312,-18.62088 34.6201,-32.207623 57.9731,-35.897619 19.5714,-4.408067 39.5199,-4.455217 59.3132,-2.055584 19.7933,2.399634 39.4313,7.246051 58.3817,12.62512 1.5096,19.173023 5.7329,49.670133 -0.3128,65.430463 -16.657,-0.26915 -39.1842,5.88414 -31.6548,-19 -1.6482,-23.82053 -21.3427,-36.82701 -42.4938,-38.50012 -21.1511,-1.67312 -43.7588,7.98715 -51.2335,29.50012 -14.1826,27.84603 4.1914,57.9202 30,70 19.9649,13.49423 43.5526,24.76438 64.3552,38.77589 20.8027,14.01152 38.8202,30.7644 47.6448,55.22411 8.8577,24.73667 4.7321,50.90082 -7.7951,72.42372 -12.5272,21.5229 -33.456,38.40456 -58.2049,44.57628 -25.5516,8.81359 -53.4285,7.61391 -80,6 z m -502.75004,-4.49583 c -14.49634,-20.16798 -29.05177,-40.31989 -43.66494,-60.44416 -14.61316,-20.12428 -29.28404,-40.22091 -44.01126,-60.27834 -14.72723,-20.05744 -29.51079,-40.07567 -44.34931,-60.04313 -14.83851,-19.96747 -29.73199,-39.88417 -44.67904,-59.73854 0.085,67.94114 -2.59454,136.28635 1.45455,204 2.37361,21.33297 41.95968,2.60802 32,29.23572 -5.68484,10.73778 -30.97918,1.86608 -44.5036,4.76428 -19.16547,0 -38.33093,0 -57.4964,0 -12.17745,-31.99983 32.34126,-11.09351 34,-38 0.9227,-19.95147 1.50618,-39.93087 1.84747,-59.92634 0.34129,-19.99547 0.44039,-40.00702 0.39435,-60.02281 -0.0461,-20.01579 -0.23725,-40.03582 -0.47657,-60.04824 -0.23931,-20.01242 -0.52674,-40.01724 -0.76525,-60.00261 10.20299,-32.81396 -63.88039,-18.66932 -25.06566,-45.731739 23.84724,1.313312 50.33669,-2.975416 73.00768,1.582846 14.56739,20.242453 29.23861,40.419643 43.95531,60.568633 14.71671,20.14899 29.47889,40.26977 44.2282,60.39943 14.74931,20.12965 29.48574,40.26817 44.15094,60.45262 14.66519,20.18445 29.25916,40.41483 43.72353,60.72821 0.83332,-70.07002 2.82145,-140.78692 -1,-210.54546 -10.39601,-9.71246 -51.54014,-21.85965 -24.14642,-33.454539 30.9003,0.49167 61.92696,-1.293841 92.71785,1.42857 12.44918,33.638319 -51.66003,12.223739 -34.00692,55.415179 -0.24121,22.59532 -0.38824,45.19117 -0.47532,67.78734 -0.0871,22.59618 -0.11419,45.19267 -0.11559,67.7893 -10e-4,22.59662 0.0229,45.19336 0.0387,67.79003 0.0158,22.59667 0.0231,45.19326 -0.0123,67.78958 -15.5356,-0.53883 -31.43722,1.23098 -46.75,-1.49583 z M 132.38637,405.67046 c -2.30195,-19.38098 35.20853,-6.01885 35.99999,-32 10.6943,-23.20341 21.35806,-46.41833 32.00899,-69.63818 10.65093,-23.21985 21.28905,-46.44461 31.93207,-69.6677 10.64303,-23.22308 21.29096,-46.44449 31.96154,-69.65761 10.67058,-23.21313 21.3638,-46.41797 32.0974,-69.607939 14.00633,-14.784664 31.06301,-3.589508 34,14.571429 9.47589,20.50998 18.90267,41.04032 28.32168,61.57386 9.41901,20.53355 18.83027,41.07031 28.27512,61.59315 9.44485,20.52284 18.92331,41.03177 28.47672,61.50964 9.55341,20.47788 19.18179,40.9247 28.92648,61.32335 8.94847,18.57494 17.97684,44.59976 44,40 6.95283,23.21831 -11.22619,21.30477 -28.79645,20 -17.53393,0 -35.06785,0 -52.60178,0 -17.53392,0 -35.06785,0 -52.60177,0 -13.61053,-31.56953 34.86478,-11.43099 30,-37 -7.96326,-20.44937 -16.90067,-40.88219 -27.17188,-60.1481 -22.69835,-0.30159 -45.74783,-1.00478 -68.70149,-1.13746 -22.95367,-0.13269 -45.81153,0.30513 -68.12663,2.28556 -10.15614,23.48394 -27.65652,46.62921 -26.57143,72.57143 21.87983,-9.59529 37.44237,23.82765 14.09999,23.42857 -25.17618,0 -50.35237,0 -75.52855,0 0,-3.33333 0,-6.66667 0,-10 z m 213.99999,-112 c -9.85189,-19.92913 -19.12968,-40.10864 -28.4649,-60.26037 -9.33521,-20.15173 -18.72785,-40.27567 -28.8094,-60.09366 -8.0471,4.62582 -13.12108,26.97607 -19.7257,38.35403 -12.60829,27.21944 -25.26661,54.41816 -37,82 18.19346,1.65846 37.60009,2.48757 57.00505,2.48751 19.40495,-6e-5 38.80822,-0.82929 56.99495,-2.48751 z m 582,112 c -1.2337,-21.88184 52.78642,-3.94103 39.43532,-42.15625 0.24195,-20.8192 0.38925,-41.63897 0.47635,-62.45909 0.0871,-20.82012 0.11398,-41.64058 0.1151,-62.46118 10e-4,-20.8206 -0.0235,-41.64134 -0.0394,-62.46199 -0.0159,-20.82065 -0.0232,-41.64122 0.0127,-62.46149 -26.37832,2.29828 -55.93232,-6.20371 -79.85787,6.80589 -1.50256,28.04473 -4.68344,46.64651 -35.23142,39.19411 -7.11643,-15.54116 -0.90104,-46.83584 -1.48214,-66.571429 25.02495,-0.70812 50.0617,-1.115779 75.10565,-1.337873 25.04395,-0.222093 50.0951,-0.258621 75.14885,-0.224476 25.0538,0.03414 50.1101,0.138961 75.1645,0.199555 25.0544,0.06059 50.1067,0.07697 75.1525,-0.06578 0,22.666673 0,45.333333 0,68.000003 -11.5708,-3.09584 -38.1128,8.77945 -32,-12 6.3763,-41.53883 -40.2943,-34.36261 -67.7072,-34 -23.6002,-7.69784 -14.9407,16.86409 -16.2929,31.41431 0.1745,19.75739 0.1467,39.53024 0.067,59.30763 -0.08,19.77738 -0.2116,39.5593 -0.2452,59.33481 -0.034,19.77551 0.031,39.54462 0.3439,59.29638 0.3129,19.75176 0.8743,39.48618 1.8343,59.19232 6.3204,22.76756 50.3983,0.0574 38.0001,33.45455 -24.6667,0 -49.3333,0 -74,0 -24.66671,0 -49.3334,0 -74.00007,0 0,-3.33333 0,-6.66667 0,-10 z";

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
var textPointsId = [];

if (textMesh){

    var myText = svgToPoints(svgString);
    points = myText.points;
    citySet = new Set(range(0, points.length));
    var scale = 0.5

    // scale points to [0,1] + scale
    var maxX = Math.max.apply(Math, points.map(function(p){return p.x}));
    var minX = Math.min.apply(Math, points.map(function(p){return p.x}));
    var maxY = Math.max.apply(Math, points.map(function(p){return p.y}));
    var minY = Math.min.apply(Math, points.map(function(p){return p.y}));
    points = points.map(function(p){
        textPointsId.push(p.id);
        return {
            id:p.id, 
            x: 0.4*(p.x-minX)/(maxX-minX)+0.25, 
            y: 0.4*(p.y-minY)/(maxY-minY)+0.25}
    });

    //add random points
    var nbPoints = points.length;
    for(var i=0; i<nbRandomPoints; ++i) {
        points.push({
            id : nbPoints,
            // x: random() * w + xInit,
            // y: random() * h + yInit
            x:random(),
            y:random()
        });
        nbPoints++;
    }

} else {
    //add random points
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
var nextEdges = new Map();
var edges = [];
//var permutations = [[0,1], [1,0], [0,2], [2,0], [1,2], [2,1]];
var permutations = [[0,1], [0,2], [1,2]];
var nbEdges = 0;
cells.forEach(function(cell){
   
    for (var i = 0; i < 3; ++i){  // for each point.id listed in current cell
        var pt = cell[i];
        var nexts = [];

        //console.log("-------- id:" + pt);

        for (var j = 1; j <= 2; ++j){ 

            var ptj = cell[( i + j ) % 3];
            //console.log("other: " + ptj);
            var temp = []; // equivalent of nexts, but for point ptj
            var newEdge = undefined;

            if ( nextEdges.has(points[pt]) ){ // if key already exists, add the corresponding edges
                //console.log("Trouvé 1");
                var t = nextEdges.get(points[pt]);
                var temppoints = t.map(function(e){
                    return [e.pt1, e.pt2];
                })

                temppoints = temppoints.reduce(function(a, b){
                     return a.concat(b);
                });

                //console.log(temppoints);

                if (temppoints.indexOf(points[ptj]) == -1){
                    //console.log("Trouvé 2");
                    newEdge = createEdge(points[pt], points[ptj]);
                    nexts.push(newEdge);
                    edges.push(newEdge);
                    //console.log(edges.length);
                }
            }
            else {
                //console.log("Trouvé 3");
                newEdge = createEdge(points[pt], points[ptj]);
                nexts.push(newEdge);
                edges.push(newEdge);
                //console.log(edges.length);
            }
            if (newEdge != undefined){
                temp.push(newEdge);
                addToNextEdges(points[ptj], cell, temp, nextEdges); // add also the edge to the edge's other point's nextEdges
            }
        }

        if (textMesh && textPointsId.indexOf(pt) != -1 && pt < (textPointsId.length - 1)) // add the textEdges to nextEdges map
        {
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

    // // edges
    // context.strokeStyle = "#000";
    // for(var i=0; i<edges.length; ++i) {
        
    //     var edge = edges[i];
    //     if (edge.pheromon != 0){
    //         context.lineWidth = Math.min(0.00001 * edge.pheromon, 0.01);
    //     } else {
    //         context.lineWidth = 0.00001;
    //     }
    //     context.beginPath();
    //     context.moveTo(points[edge.pt1.id].x, points[edge.pt1.id].y);
    //     context.lineTo(points[edge.pt2.id].x, points[edge.pt2.id].y);
    //     context.stroke();
    // }

    // // vertices
    // for(var i=0; i<points.length; ++i) {
    //     context.beginPath()
    //     var point = points[i];
    //     if (citySet.has(point.id)) {
    //         context.fillStyle = "#0101DF";
    //         context.arc(point.x, point.y, 0.006, 0, 2*Math.PI);
    //     }
    //     else {
    //         context.fillStyle = "#000";
    //         context.arc(points[i].x, points[i].y, 0.003, 0, 2*Math.PI);
    //     }
    //     context.closePath();
    //     context.fill();
    // }

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
                        // increased dropped pheromons for textEdges
                        if ((textPointsId.indexOf(a.id)) != -1 && (textPointsId.indexOf(b.id) != -1) && (Math.abs(a.id - b.id) == 1))
                        {
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