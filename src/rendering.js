'use strict'

var antFunction = require('./ant.js');

var random = Math.random;

var RANDOMMVT = 0.003;
var ANTSIZE = 0.002;

module.exports = function(container, initVar){
	
	if(!container)
		throw new TypeError('Missing container');

	var edges = initVar.edges;
	var population = initVar.population;
	var pointsInfos = initVar.pointsInfos;
	var nbAnts = population.length;

	var canvasList = document.getElementsByTagName("canvas");
	
	if (canvasList.length === 0){
		var canvas = document.createElement("canvas");
		var rect = container.getBoundingClientRect();
		canvas.width = rect.width;
		canvas.height = rect.height;
		canvas.style.backgroundColor = "rgba(250, 250, 250, 0)"; 
		container.appendChild(canvas);
	}
	else{
		var canvas = canvasList[0];
		console.log('CANVAS');
	}
	
	var context = canvas.getContext("2d");
	context.clearRect ( 0 , 0 , canvas.width, canvas.height );
	

	function tick() {
		var w = canvas.width;
		var h = canvas.height;
		var mouse = [lastMouseMoveEvent.clientX/w, lastMouseMoveEvent.clientY/h];
		context.setTransform(w, 0, 0, h, 0, 0);
		context.fillStyle = "rgba(250, 250, 250, 0.4)";
		context.fillRect(0,0,w,h);

		// edges
		// context.strokeStyle = "#000";
		// for(var i=0; i < edges.length; ++i) {
		//     context.lineWidth = 0.0001;
		//     var edge = edges[i];
		//     // if (edge.pheromon != 0){
		//     //     context.lineWidth = Math.min(0.00001 * edge.pheromon, 0.01);
		//     // } else {
		//     //     context.lineWidth = 0.00001;
		//     // }
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
		population.forEach(function(ant){
			ant.transit();
		});

		// for (i = 0; i < nbAnts; i++) {
		//     population[i].transit();
		// }

		// pheromon evaporation
		edges.forEach(function(edge){
			if(edge.pheromon > 0){
				edge.pheromon -= 0.0001;
			}
		});

		// for (i = 0; i < edges.length; i++) {
			

		// ants
		population.forEach(function(ant){
			context.beginPath()
			var x = ant.x + RANDOMMVT*random();
			var y = ant.y + RANDOMMVT*random();

			context.fillStyle = "black"
			context.fillRect(x, y, ANTSIZE, ANTSIZE);
			context.closePath();
			context.fill();
		})
		// for(var i=0; i<population.length; ++i) {
			
		// }

	};
	
	var lastMouseMoveEvent = {
		clientX: 0,
		clientY: 0
	};
	
	container.addEventListener('mousemove', function(e){
		lastMouseMoveEvent = e;
	});
	
	var paused = false;
	
	function togglePlayPause(){
		paused = !paused;
		if(!paused)
			animate();
	}
	
	// container.addEventListener('click', togglePlayPause);

	function animate(){
		tick();
		
		if(!paused)
			requestAnimationFrame(animate);
	}
	animate();


	function setAntCount(opts){

		var previousCount = population.length;

		if (opts.nbAnts > previousCount){
			var Ant = antFunction(container, pointsInfos, opts);

			for (var i = 0; i < options.nbAnts - previousCount; i++) {
				var newAnt = new Ant(Ant.generateRandStartPoint());
				newAnt.setDirection();
				population.push = newAnt;
			}
		}
		else{
			population = population.slice(0, opts.nbAnts);
			nbAnts = population.length;
			console.log('Nb Ants :', population.length);
		}
		
	}

	function modifyAnts(opts){
		setAntCount(opts);

		population.forEach(function(ant){
			ant.velocity = opts.velocity;
			ant.weight = opts.weight;
			ant.repSize = opts.repSize;
			ant.repSpeed = opts.repSpeed;
		});
	}
	
	return {
		togglePlayPause: togglePlayPause,
		// should be a getter/setter, but IE8
		getAntCount: function(){
			return population.length;
		},
		setAntCount: setAntCount,
		modifyAnts: modifyAnts
	}
}
