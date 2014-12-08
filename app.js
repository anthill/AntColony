(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/index.js":[function(require,module,exports){
'use strict';

var initRendering = require('./src/rendering.js');
var initializePoints = require('./src/initializePoints.js');
var createEdges = require('./src/createEdges.js');
var initAnts = require('./src/initializeAnts');

module.exports = function init(containerElement, options){

	var render, pointsInfos, edges, population, initVar;


	function _init(containerElement, options){
		pointsInfos = initializePoints(options.nbStart, options.nbRand);
		edges = createEdges(pointsInfos.points);
		population = initAnts(containerElement, pointsInfos, options);
		initVar = {
			pointsInfos: pointsInfos,
			edges: edges,
			population: population
		};
		render = initRendering(containerElement, initVar);
	}

	_init(containerElement, options);

	return {
		togglePlayPause: function(){ render.togglePlayPause() },
		changeOptions: function(opts){
			if (opts.nbStart === undefined && opts.nbRand === undefined){
				// modify population
				render.modifyAnts(opts);
			}
			else{
				// reset previous animation
				render.reset();

				// reset elements
				_init(containerElement, opts);
				
			}
		},
		reset: function(opts){
			render.reset();

				// reset elements
			_init(containerElement, opts);
		}
	};
};
},{"./src/createEdges.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/createEdges.js","./src/initializeAnts":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/initializeAnts.js","./src/initializePoints.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/initializePoints.js","./src/rendering.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/rendering.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js":[function(require,module,exports){
"use strict"

//High level idea:
// 1. Use Clarkson's incremental construction to find convex hull
// 2. Point location in triangulation by jump and walk

module.exports = incrementalConvexHull

var orient = require("robust-orientation")
var compareCell = require("simplicial-complex").compareCells

function compareInt(a, b) {
  return a - b
}

function Simplex(vertices, adjacent, boundary) {
  this.vertices = vertices
  this.adjacent = adjacent
  this.boundary = boundary
  this.lastVisited = -1
}

Simplex.prototype.flip = function() {
  var t = this.vertices[0]
  this.vertices[0] = this.vertices[1]
  this.vertices[1] = t
  var u = this.adjacent[0]
  this.adjacent[0] = this.adjacent[1]
  this.adjacent[1] = u
}

function GlueFacet(vertices, cell, index) {
  this.vertices = vertices
  this.cell = cell
  this.index = index
}

function compareGlue(a, b) {
  return compareCell(a.vertices, b.vertices)
}

function bakeOrient(d) {
  var code = ["function orient(){var tuple=this.tuple;return test("]
  for(var i=0; i<=d; ++i) {
    if(i > 0) {
      code.push(",")
    }
    code.push("tuple[", i, "]")
  }
  code.push(")}return orient")
  var proc = new Function("test", code.join(""))
  var test = orient[d+1]
  if(!test) {
    test = orient
  }
  return proc(test)
}

var BAKED = []

function Triangulation(dimension, vertices, simplices) {
  this.dimension = dimension
  this.vertices = vertices
  this.simplices = simplices
  this.interior = simplices.filter(function(c) {
    return !c.boundary
  })

  this.tuple = new Array(dimension+1)
  for(var i=0; i<=dimension; ++i) {
    this.tuple[i] = this.vertices[i]
  }

  var o = BAKED[dimension]
  if(!o) {
    o = BAKED[dimension] = bakeOrient(dimension)
  }
  this.orient = o
}

var proto = Triangulation.prototype

//Degenerate situation where we are on boundary, but coplanar to face
proto.handleBoundaryDegeneracy = function(cell, point) {
  var d = this.dimension
  var n = this.vertices.length - 1
  var tuple = this.tuple
  var verts = this.vertices

  //Dumb solution: Just do dfs from boundary cell until we find any peak, or terminate
  var toVisit = [ cell ]
  cell.lastVisited = -n
  while(toVisit.length > 0) {
    cell = toVisit.pop()
    var cellVerts = cell.vertices
    var cellAdj = cell.adjacent
    for(var i=0; i<=d; ++i) {
      var neighbor = cellAdj[i]
      if(!neighbor.boundary || neighbor.lastVisited <= -n) {
        continue
      }
      var nv = neighbor.vertices
      for(var j=0; j<=d; ++j) {
        var vv = nv[j]
        if(vv < 0) {
          tuple[j] = point
        } else {
          tuple[j] = verts[vv]
        }
      }
      var o = this.orient()
      if(o > 0) {
        return neighbor
      }
      neighbor.lastVisited = -n
      if(o === 0) {
        toVisit.push(neighbor)
      }
    }
  }
  return null
}

proto.walk = function(point, random) {
  //Alias local properties
  var n = this.vertices.length - 1
  var d = this.dimension
  var verts = this.vertices
  var tuple = this.tuple

  //Compute initial jump cell
  var initIndex = random ? (this.interior.length * Math.random())|0 : (this.interior.length-1)
  var cell = this.interior[ initIndex ]

  //Start walking
outerLoop:
  while(!cell.boundary) {
    var cellVerts = cell.vertices
    var cellAdj = cell.adjacent

    for(var i=0; i<=d; ++i) {
      tuple[i] = verts[cellVerts[i]]
    }
    cell.lastVisited = n

    //Find farthest adjacent cell
    for(var i=0; i<=d; ++i) {
      var neighbor = cellAdj[i]
      if(neighbor.lastVisited >= n) {
        continue
      }
      var prev = tuple[i]
      tuple[i] = point
      var o = this.orient()
      tuple[i] = prev
      if(o < 0) {
        cell = neighbor
        continue outerLoop
      } else {
        if(!neighbor.boundary) {
          neighbor.lastVisited = n
        } else {
          neighbor.lastVisited = -n
        }
      }
    }
    return
  }

  return cell
}

proto.addPeaks = function(point, cell) {
  var n = this.vertices.length - 1
  var d = this.dimension
  var verts = this.vertices
  var tuple = this.tuple
  var interior = this.interior
  var simplices = this.simplices

  //Walking finished at boundary, time to add peaks
  var tovisit = [ cell ]

  //Stretch initial boundary cell into a peak
  cell.lastVisited = n
  cell.vertices[cell.vertices.indexOf(-1)] = n
  cell.boundary = false
  interior.push(cell)

  //Record a list of all new boundaries created by added peaks so we can glue them together when we are all done
  var glueFacets = []

  //Do a traversal of the boundary walking outward from starting peak
  while(tovisit.length > 0) {
    //Pop off peak and walk over adjacent cells
    var cell = tovisit.pop()
    var cellVerts = cell.vertices
    var cellAdj = cell.adjacent
    var indexOfN = cellVerts.indexOf(n)
    if(indexOfN < 0) {
      continue
    }

    for(var i=0; i<=d; ++i) {
      if(i === indexOfN) {
        continue
      }

      //For each boundary neighbor of the cell
      var neighbor = cellAdj[i]
      if(!neighbor.boundary || neighbor.lastVisited >= n) {
        continue
      }

      var nv = neighbor.vertices

      //Test if neighbor is a peak
      if(neighbor.lastVisited !== -n) {      
        //Compute orientation of p relative to each boundary peak
        var indexOfNeg1 = 0
        for(var j=0; j<=d; ++j) {
          if(nv[j] < 0) {
            indexOfNeg1 = j
            tuple[j] = point
          } else {
            tuple[j] = verts[nv[j]]
          }
        }
        var o = this.orient()

        //Test if neighbor cell is also a peak
        if(o > 0) {
          nv[indexOfNeg1] = n
          neighbor.boundary = false
          interior.push(neighbor)
          tovisit.push(neighbor)
          neighbor.lastVisited = n
          continue
        } else {
          neighbor.lastVisited = -n
        }
      }

      var na = neighbor.adjacent

      //Otherwise, replace neighbor with new face
      var vverts = cellVerts.slice()
      var vadj = cellAdj.slice()
      var ncell = new Simplex(vverts, vadj, true)
      simplices.push(ncell)

      //Connect to neighbor
      var opposite = na.indexOf(cell)
      if(opposite < 0) {
        continue
      }
      na[opposite] = ncell
      vadj[indexOfN] = neighbor

      //Connect to cell
      vverts[i] = -1
      vadj[i] = cell
      cellAdj[i] = ncell

      //Flip facet
      ncell.flip()

      //Add to glue list
      for(var j=0; j<=d; ++j) {
        var uu = vverts[j]
        if(uu < 0 || uu === n) {
          continue
        }
        var nface = new Array(d-1)
        var nptr = 0
        for(var k=0; k<=d; ++k) {
          var vv = vverts[k]
          if(vv < 0 || k === j) {
            continue
          }
          nface[nptr++] = vv
        }
        glueFacets.push(new GlueFacet(nface, ncell, j))
      }
    }
  }

  //Glue boundary facets together
  glueFacets.sort(compareGlue)

  for(var i=0; i+1<glueFacets.length; i+=2) {
    var a = glueFacets[i]
    var b = glueFacets[i+1]
    var ai = a.index
    var bi = b.index
    if(ai < 0 || bi < 0) {
      continue
    }
    a.cell.adjacent[a.index] = b.cell
    b.cell.adjacent[b.index] = a.cell
  }
}

proto.insert = function(point, random) {
  //Add point
  var verts = this.vertices
  verts.push(point)

  var cell = this.walk(point, random)
  if(!cell) {
    return
  }

  //Alias local properties
  var d = this.dimension
  var tuple = this.tuple

  //Degenerate case: If point is coplanar to cell, then walk until we find a non-degenerate boundary
  for(var i=0; i<=d; ++i) {
    var vv = cell.vertices[i]
    if(vv < 0) {
      tuple[i] = point
    } else {
      tuple[i] = verts[vv]
    }
  }
  var o = this.orient(tuple)
  if(o < 0) {
    return
  } else if(o === 0) {
    cell = this.handleBoundaryDegeneracy(cell, point)
    if(!cell) {
      return
    }
  }

  //Add peaks
  this.addPeaks(point, cell)
}

//Extract all boundary cells
proto.boundary = function() {
  var d = this.dimension
  var boundary = []
  var cells = this.simplices
  var nc = cells.length
  for(var i=0; i<nc; ++i) {
    var c = cells[i]
    if(c.boundary) {
      var bcell = new Array(d)
      var cv = c.vertices
      var ptr = 0
      var parity = 0
      for(var j=0; j<=d; ++j) {
        if(cv[j] >= 0) {
          bcell[ptr++] = cv[j]
        } else {
          parity = j&1
        }
      }
      if(parity === (d&1)) {
        var t = bcell[0]
        bcell[0] = bcell[1]
        bcell[1] = t
      }
      boundary.push(bcell)
    }
  }
  return boundary
}

function incrementalConvexHull(points, randomSearch) {
  var n = points.length
  if(n === 0) {
    throw new Error("Must have at least d+1 points")
  }
  var d = points[0].length
  if(n <= d) {
    throw new Error("Must input at least d+1 points")
  }

  //FIXME: This could be degenerate, but need to select d+1 non-coplanar points to bootstrap process
  var initialSimplex = points.slice(0, d+1)

  //Make sure initial simplex is positively oriented
  var o = orient.apply(void 0, initialSimplex)
  if(o === 0) {
    throw new Error("Input not in general position")
  }
  var initialCoords = new Array(d+1)
  for(var i=0; i<=d; ++i) {
    initialCoords[i] = i
  }
  if(o < 0) {
    initialCoords[0] = 1
    initialCoords[1] = 0
  }

  //Create initial topological index, glue pointers together (kind of messy)
  var initialCell = new Simplex(initialCoords, new Array(d+1), false)
  var boundary = initialCell.adjacent
  var list = new Array(d+2)
  for(var i=0; i<=d; ++i) {
    var verts = initialCoords.slice()
    for(var j=0; j<=d; ++j) {
      if(j === i) {
        verts[j] = -1
      }
    }
    var t = verts[0]
    verts[0] = verts[1]
    verts[1] = t
    var cell = new Simplex(verts, new Array(d+1), true)
    boundary[i] = cell
    list[i] = cell
  }
  list[d+1] = initialCell
  for(var i=0; i<=d; ++i) {
    var verts = boundary[i].vertices
    var adj = boundary[i].adjacent
    for(var j=0; j<=d; ++j) {
      var v = verts[j]
      if(v < 0) {
        adj[j] = initialCell
        continue
      }
      for(var k=0; k<=d; ++k) {
        if(boundary[k].vertices.indexOf(v) < 0) {
          adj[j] = boundary[k]
        }
      }
    }
  }

  //Initialize triangles
  var triangles = new Triangulation(d, initialSimplex, list)

  //Insert remaining points
  var useRandom = !!randomSearch
  for(var i=d+1; i<n; ++i) {
    triangles.insert(points[i], useRandom)
  }
  
  //Extract boundary cells
  return triangles.boundary()
}
},{"robust-orientation":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js","simplicial-complex":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js":[function(require,module,exports){
"use strict"

module.exports = fastTwoSum

function fastTwoSum(a, b, result) {
	var x = a + b
	var bv = x - a
	var av = x - bv
	var br = b - bv
	var ar = a - av
	if(result) {
		result[0] = ar + br
		result[1] = x
		return result
	}
	return [ar+br, x]
}
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js":[function(require,module,exports){
"use strict"

var twoProduct = require("two-product")
var twoSum = require("two-sum")

module.exports = scaleLinearExpansion

function scaleLinearExpansion(e, scale) {
  var n = e.length
  if(n === 1) {
    var ts = twoProduct(e[0], scale)
    if(ts[0]) {
      return ts
    }
    return [ ts[1] ]
  }
  var g = new Array(2 * n)
  var q = [0.1, 0.1]
  var t = [0.1, 0.1]
  var count = 0
  twoProduct(e[0], scale, q)
  if(q[0]) {
    g[count++] = q[0]
  }
  for(var i=1; i<n; ++i) {
    twoProduct(e[i], scale, t)
    var pq = q[1]
    twoSum(pq, t[0], q)
    if(q[0]) {
      g[count++] = q[0]
    }
    var a = t[1]
    var b = q[1]
    var x = a + b
    var bv = x - a
    var y = b - bv
    q[1] = x
    if(y) {
      g[count++] = y
    }
  }
  if(q[1]) {
    g[count++] = q[1]
  }
  if(count === 0) {
    g[count++] = 0.0
  }
  g.length = count
  return g
}
},{"two-product":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js","two-sum":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js":[function(require,module,exports){
"use strict"

module.exports = robustSubtract

//Easy case: Add two scalars
function scalarScalar(a, b) {
  var x = a + b
  var bv = x - a
  var av = x - bv
  var br = b - bv
  var ar = a - av
  var y = ar + br
  if(y) {
    return [y, x]
  }
  return [x]
}

function robustSubtract(e, f) {
  var ne = e.length|0
  var nf = f.length|0
  if(ne === 1 && nf === 1) {
    return scalarScalar(e[0], -f[0])
  }
  var n = ne + nf
  var g = new Array(n)
  var count = 0
  var eptr = 0
  var fptr = 0
  var abs = Math.abs
  var ei = e[eptr]
  var ea = abs(ei)
  var fi = -f[fptr]
  var fa = abs(fi)
  var a, b
  if(ea < fa) {
    b = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    b = fi
    fptr += 1
    if(fptr < nf) {
      fi = -f[fptr]
      fa = abs(fi)
    }
  }
  if((eptr < ne && ea < fa) || (fptr >= nf)) {
    a = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    a = fi
    fptr += 1
    if(fptr < nf) {
      fi = -f[fptr]
      fa = abs(fi)
    }
  }
  var x = a + b
  var bv = x - a
  var y = b - bv
  var q0 = y
  var q1 = x
  var _x, _bv, _av, _br, _ar
  while(eptr < ne && fptr < nf) {
    if(ea < fa) {
      a = ei
      eptr += 1
      if(eptr < ne) {
        ei = e[eptr]
        ea = abs(ei)
      }
    } else {
      a = fi
      fptr += 1
      if(fptr < nf) {
        fi = -f[fptr]
        fa = abs(fi)
      }
    }
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
  }
  while(eptr < ne) {
    a = ei
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
    }
  }
  while(fptr < nf) {
    a = fi
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    } 
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    fptr += 1
    if(fptr < nf) {
      fi = -f[fptr]
    }
  }
  if(q0) {
    g[count++] = q0
  }
  if(q1) {
    g[count++] = q1
  }
  if(!count) {
    g[count++] = 0.0  
  }
  g.length = count
  return g
}
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js":[function(require,module,exports){
"use strict"

module.exports = linearExpansionSum

//Easy case: Add two scalars
function scalarScalar(a, b) {
  var x = a + b
  var bv = x - a
  var av = x - bv
  var br = b - bv
  var ar = a - av
  var y = ar + br
  if(y) {
    return [y, x]
  }
  return [x]
}

function linearExpansionSum(e, f) {
  var ne = e.length|0
  var nf = f.length|0
  if(ne === 1 && nf === 1) {
    return scalarScalar(e[0], f[0])
  }
  var n = ne + nf
  var g = new Array(n)
  var count = 0
  var eptr = 0
  var fptr = 0
  var abs = Math.abs
  var ei = e[eptr]
  var ea = abs(ei)
  var fi = f[fptr]
  var fa = abs(fi)
  var a, b
  if(ea < fa) {
    b = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    b = fi
    fptr += 1
    if(fptr < nf) {
      fi = f[fptr]
      fa = abs(fi)
    }
  }
  if((eptr < ne && ea < fa) || (fptr >= nf)) {
    a = ei
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
      ea = abs(ei)
    }
  } else {
    a = fi
    fptr += 1
    if(fptr < nf) {
      fi = f[fptr]
      fa = abs(fi)
    }
  }
  var x = a + b
  var bv = x - a
  var y = b - bv
  var q0 = y
  var q1 = x
  var _x, _bv, _av, _br, _ar
  while(eptr < ne && fptr < nf) {
    if(ea < fa) {
      a = ei
      eptr += 1
      if(eptr < ne) {
        ei = e[eptr]
        ea = abs(ei)
      }
    } else {
      a = fi
      fptr += 1
      if(fptr < nf) {
        fi = f[fptr]
        fa = abs(fi)
      }
    }
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
  }
  while(eptr < ne) {
    a = ei
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    }
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    eptr += 1
    if(eptr < ne) {
      ei = e[eptr]
    }
  }
  while(fptr < nf) {
    a = fi
    b = q0
    x = a + b
    bv = x - a
    y = b - bv
    if(y) {
      g[count++] = y
    } 
    _x = q1 + x
    _bv = _x - q1
    _av = _x - _bv
    _br = x - _bv
    _ar = q1 - _av
    q0 = _ar + _br
    q1 = _x
    fptr += 1
    if(fptr < nf) {
      fi = f[fptr]
    }
  }
  if(q0) {
    g[count++] = q0
  }
  if(q1) {
    g[count++] = q1
  }
  if(!count) {
    g[count++] = 0.0  
  }
  g.length = count
  return g
}
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js":[function(require,module,exports){
"use strict"

module.exports = twoProduct

var SPLITTER = +(Math.pow(2, 27) + 1.0)

function twoProduct(a, b, result) {
  var x = a * b

  var c = SPLITTER * a
  var abig = c - a
  var ahi = c - abig
  var alo = a - ahi

  var d = SPLITTER * b
  var bbig = d - b
  var bhi = d - bbig
  var blo = b - bhi

  var err1 = x - (ahi * bhi)
  var err2 = err1 - (alo * bhi)
  var err3 = err2 - (ahi * blo)

  var y = alo * blo - err3

  if(result) {
    result[0] = y
    result[1] = x
    return result
  }

  return [ y, x ]
}
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js":[function(require,module,exports){
"use strict"

var twoProduct = require("two-product")
var robustSum = require("robust-sum")
var robustScale = require("robust-scale")
var robustSubtract = require("robust-subtract")

var NUM_EXPAND = 5

var EPSILON     = 1.1102230246251565e-16
var ERRBOUND3   = (3.0 + 16.0 * EPSILON) * EPSILON
var ERRBOUND4   = (7.0 + 56.0 * EPSILON) * EPSILON

function cofactor(m, c) {
  var result = new Array(m.length-1)
  for(var i=1; i<m.length; ++i) {
    var r = result[i-1] = new Array(m.length-1)
    for(var j=0,k=0; j<m.length; ++j) {
      if(j === c) {
        continue
      }
      r[k++] = m[i][j]
    }
  }
  return result
}

function matrix(n) {
  var result = new Array(n)
  for(var i=0; i<n; ++i) {
    result[i] = new Array(n)
    for(var j=0; j<n; ++j) {
      result[i][j] = ["m", j, "[", (n-i-1), "]"].join("")
    }
  }
  return result
}

function sign(n) {
  if(n & 1) {
    return "-"
  }
  return ""
}

function generateSum(expr) {
  if(expr.length === 1) {
    return expr[0]
  } else if(expr.length === 2) {
    return ["sum(", expr[0], ",", expr[1], ")"].join("")
  } else {
    var m = expr.length>>1
    return ["sum(", generateSum(expr.slice(0, m)), ",", generateSum(expr.slice(m)), ")"].join("")
  }
}

function determinant(m) {
  if(m.length === 2) {
    return [["sum(prod(", m[0][0], ",", m[1][1], "),prod(-", m[0][1], ",", m[1][0], "))"].join("")]
  } else {
    var expr = []
    for(var i=0; i<m.length; ++i) {
      expr.push(["scale(", generateSum(determinant(cofactor(m, i))), ",", sign(i), m[0][i], ")"].join(""))
    }
    return expr
  }
}

function orientation(n) {
  var pos = []
  var neg = []
  var m = matrix(n)
  var args = []
  for(var i=0; i<n; ++i) {
    if((i&1)===0) {
      pos.push.apply(pos, determinant(cofactor(m, i)))
    } else {
      neg.push.apply(neg, determinant(cofactor(m, i)))
    }
    args.push("m" + i)
  }
  var posExpr = generateSum(pos)
  var negExpr = generateSum(neg)
  var funcName = "orientation" + n + "Exact"
  var code = ["function ", funcName, "(", args.join(), "){var p=", posExpr, ",n=", negExpr, ",d=sub(p,n);\
return d[d.length-1];};return ", funcName].join("")
  var proc = new Function("sum", "prod", "scale", "sub", code)
  return proc(robustSum, twoProduct, robustScale, robustSubtract)
}

var orientation3Exact = orientation(3)
var orientation4Exact = orientation(4)

var CACHED = [
  function orientation0() { return 0 },
  function orientation1() { return 0 },
  function orientation2(a, b) { 
    return b[0] - a[0]
  },
  function orientation3(a, b, c) {
    var l = (a[1] - c[1]) * (b[0] - c[0])
    var r = (a[0] - c[0]) * (b[1] - c[1])
    var det = l - r
    var s
    if(l > 0) {
      if(r <= 0) {
        return det
      } else {
        s = l + r
      }
    } else if(l < 0) {
      if(r >= 0) {
        return det
      } else {
        s = -(l + r)
      }
    } else {
      return det
    }
    var tol = ERRBOUND3 * s
    if(det >= tol || det <= -tol) {
      return det
    }
    return orientation3Exact(a, b, c)
  },
  function orientation4(a,b,c,d) {
    var adx = a[0] - d[0]
    var bdx = b[0] - d[0]
    var cdx = c[0] - d[0]
    var ady = a[1] - d[1]
    var bdy = b[1] - d[1]
    var cdy = c[1] - d[1]
    var adz = a[2] - d[2]
    var bdz = b[2] - d[2]
    var cdz = c[2] - d[2]
    var bdxcdy = bdx * cdy
    var cdxbdy = cdx * bdy
    var cdxady = cdx * ady
    var adxcdy = adx * cdy
    var adxbdy = adx * bdy
    var bdxady = bdx * ady
    var det = adz * (bdxcdy - cdxbdy) 
            + bdz * (cdxady - adxcdy)
            + cdz * (adxbdy - bdxady)
    var permanent = (Math.abs(bdxcdy) + Math.abs(cdxbdy)) * Math.abs(adz)
                  + (Math.abs(cdxady) + Math.abs(adxcdy)) * Math.abs(bdz)
                  + (Math.abs(adxbdy) + Math.abs(bdxady)) * Math.abs(cdz)
    var tol = ERRBOUND4 * permanent
    if ((det > tol) || (-det > tol)) {
      return det
    }
    return orientation4Exact(a,b,c,d)
  }
]

function slowOrient(args) {
  var proc = CACHED[args.length]
  if(!proc) {
    proc = CACHED[args.length] = orientation(args.length)
  }
  return proc.apply(undefined, args)
}

function generateOrientationProc() {
  while(CACHED.length <= NUM_EXPAND) {
    CACHED.push(orientation(CACHED.length))
  }
  var args = []
  var procArgs = ["slow"]
  for(var i=0; i<=NUM_EXPAND; ++i) {
    args.push("a" + i)
    procArgs.push("o" + i)
  }
  var code = [
    "function getOrientation(", args.join(), "){switch(arguments.length){case 0:case 1:return 0;"
  ]
  for(var i=2; i<=NUM_EXPAND; ++i) {
    code.push("case ", i, ":return o", i, "(", args.slice(0, i).join(), ");")
  }
  code.push("}var s=new Array(arguments.length);for(var i=0;i<arguments.length;++i){s[i]=arguments[i]};return slow(s);}return getOrientation")
  procArgs.push(code.join(""))

  var proc = Function.apply(undefined, procArgs)
  module.exports = proc.apply(undefined, [slowOrient].concat(CACHED))
  for(var i=0; i<=NUM_EXPAND; ++i) {
    module.exports[i] = CACHED[i]
  }
}

generateOrientationProc()
},{"robust-scale":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js","robust-subtract":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js","robust-sum":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js","two-product":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js":[function(require,module,exports){
/**
 * Bit twiddling hacks for JavaScript.
 *
 * Author: Mikola Lysenko
 *
 * Ported from Stanford bit twiddling hack library:
 *    http://graphics.stanford.edu/~seander/bithacks.html
 */

"use strict"; "use restrict";

//Number of bits in an integer
var INT_BITS = 32;

//Constants
exports.INT_BITS  = INT_BITS;
exports.INT_MAX   =  0x7fffffff;
exports.INT_MIN   = -1<<(INT_BITS-1);

//Returns -1, 0, +1 depending on sign of x
exports.sign = function(v) {
  return (v > 0) - (v < 0);
}

//Computes absolute value of integer
exports.abs = function(v) {
  var mask = v >> (INT_BITS-1);
  return (v ^ mask) - mask;
}

//Computes minimum of integers x and y
exports.min = function(x, y) {
  return y ^ ((x ^ y) & -(x < y));
}

//Computes maximum of integers x and y
exports.max = function(x, y) {
  return x ^ ((x ^ y) & -(x < y));
}

//Checks if a number is a power of two
exports.isPow2 = function(v) {
  return !(v & (v-1)) && (!!v);
}

//Computes log base 2 of v
exports.log2 = function(v) {
  var r, shift;
  r =     (v > 0xFFFF) << 4; v >>>= r;
  shift = (v > 0xFF  ) << 3; v >>>= shift; r |= shift;
  shift = (v > 0xF   ) << 2; v >>>= shift; r |= shift;
  shift = (v > 0x3   ) << 1; v >>>= shift; r |= shift;
  return r | (v >> 1);
}

//Computes log base 10 of v
exports.log10 = function(v) {
  return  (v >= 1000000000) ? 9 : (v >= 100000000) ? 8 : (v >= 10000000) ? 7 :
          (v >= 1000000) ? 6 : (v >= 100000) ? 5 : (v >= 10000) ? 4 :
          (v >= 1000) ? 3 : (v >= 100) ? 2 : (v >= 10) ? 1 : 0;
}

//Counts number of bits
exports.popCount = function(v) {
  v = v - ((v >>> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
  return ((v + (v >>> 4) & 0xF0F0F0F) * 0x1010101) >>> 24;
}

//Counts number of trailing zeros
function countTrailingZeros(v) {
  var c = 32;
  v &= -v;
  if (v) c--;
  if (v & 0x0000FFFF) c -= 16;
  if (v & 0x00FF00FF) c -= 8;
  if (v & 0x0F0F0F0F) c -= 4;
  if (v & 0x33333333) c -= 2;
  if (v & 0x55555555) c -= 1;
  return c;
}
exports.countTrailingZeros = countTrailingZeros;

//Rounds to next power of 2
exports.nextPow2 = function(v) {
  v += v === 0;
  --v;
  v |= v >>> 1;
  v |= v >>> 2;
  v |= v >>> 4;
  v |= v >>> 8;
  v |= v >>> 16;
  return v + 1;
}

//Rounds down to previous power of 2
exports.prevPow2 = function(v) {
  v |= v >>> 1;
  v |= v >>> 2;
  v |= v >>> 4;
  v |= v >>> 8;
  v |= v >>> 16;
  return v - (v>>>1);
}

//Computes parity of word
exports.parity = function(v) {
  v ^= v >>> 16;
  v ^= v >>> 8;
  v ^= v >>> 4;
  v &= 0xf;
  return (0x6996 >>> v) & 1;
}

var REVERSE_TABLE = new Array(256);

(function(tab) {
  for(var i=0; i<256; ++i) {
    var v = i, r = i, s = 7;
    for (v >>>= 1; v; v >>>= 1) {
      r <<= 1;
      r |= v & 1;
      --s;
    }
    tab[i] = (r << s) & 0xff;
  }
})(REVERSE_TABLE);

//Reverse bits in a 32 bit word
exports.reverse = function(v) {
  return  (REVERSE_TABLE[ v         & 0xff] << 24) |
          (REVERSE_TABLE[(v >>> 8)  & 0xff] << 16) |
          (REVERSE_TABLE[(v >>> 16) & 0xff] << 8)  |
           REVERSE_TABLE[(v >>> 24) & 0xff];
}

//Interleave bits of 2 coordinates with 16 bits.  Useful for fast quadtree codes
exports.interleave2 = function(x, y) {
  x &= 0xFFFF;
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;

  y &= 0xFFFF;
  y = (y | (y << 8)) & 0x00FF00FF;
  y = (y | (y << 4)) & 0x0F0F0F0F;
  y = (y | (y << 2)) & 0x33333333;
  y = (y | (y << 1)) & 0x55555555;

  return x | (y << 1);
}

//Extracts the nth interleaved component
exports.deinterleave2 = function(v, n) {
  v = (v >>> n) & 0x55555555;
  v = (v | (v >>> 1))  & 0x33333333;
  v = (v | (v >>> 2))  & 0x0F0F0F0F;
  v = (v | (v >>> 4))  & 0x00FF00FF;
  v = (v | (v >>> 16)) & 0x000FFFF;
  return (v << 16) >> 16;
}


//Interleave bits of 3 coordinates, each with 10 bits.  Useful for fast octree codes
exports.interleave3 = function(x, y, z) {
  x &= 0x3FF;
  x  = (x | (x<<16)) & 4278190335;
  x  = (x | (x<<8))  & 251719695;
  x  = (x | (x<<4))  & 3272356035;
  x  = (x | (x<<2))  & 1227133513;

  y &= 0x3FF;
  y  = (y | (y<<16)) & 4278190335;
  y  = (y | (y<<8))  & 251719695;
  y  = (y | (y<<4))  & 3272356035;
  y  = (y | (y<<2))  & 1227133513;
  x |= (y << 1);
  
  z &= 0x3FF;
  z  = (z | (z<<16)) & 4278190335;
  z  = (z | (z<<8))  & 251719695;
  z  = (z | (z<<4))  & 3272356035;
  z  = (z | (z<<2))  & 1227133513;
  
  return x | (z << 2);
}

//Extracts nth interleaved component of a 3-tuple
exports.deinterleave3 = function(v, n) {
  v = (v >>> n)       & 1227133513;
  v = (v | (v>>>2))   & 3272356035;
  v = (v | (v>>>4))   & 251719695;
  v = (v | (v>>>8))   & 4278190335;
  v = (v | (v>>>16))  & 0x3FF;
  return (v<<22)>>22;
}

//Computes next combination in colexicographic order (this is mistakenly called nextPermutation on the bit twiddling hacks page)
exports.nextCombination = function(v) {
  var t = v | (v - 1);
  return (t + 1) | (((~t & -~t) - 1) >>> (countTrailingZeros(v) + 1));
}


},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js":[function(require,module,exports){
"use strict"; "use restrict";

module.exports = UnionFind;

function UnionFind(count) {
  this.roots = new Array(count);
  this.ranks = new Array(count);
  
  for(var i=0; i<count; ++i) {
    this.roots[i] = i;
    this.ranks[i] = 0;
  }
}

var proto = UnionFind.prototype

Object.defineProperty(proto, "length", {
  "get": function() {
    return this.roots.length
  }
})

proto.makeSet = function() {
  var n = this.roots.length;
  this.roots.push(n);
  this.ranks.push(0);
  return n;
}

proto.find = function(x) {
  var roots = this.roots;
  while(roots[x] !== x) {
    var y = roots[x];
    roots[x] = roots[y];
    x = y;
  }
  return x;
}

proto.link = function(x, y) {
  var xr = this.find(x)
    , yr = this.find(y);
  if(xr === yr) {
    return;
  }
  var ranks = this.ranks
    , roots = this.roots
    , xd    = ranks[xr]
    , yd    = ranks[yr];
  if(xd < yd) {
    roots[xr] = yr;
  } else if(yd < xd) {
    roots[yr] = xr;
  } else {
    roots[yr] = xr;
    ++ranks[xr];
  }
}
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js":[function(require,module,exports){
"use strict"; "use restrict";

var bits      = require("bit-twiddle")
  , UnionFind = require("union-find")

//Returns the dimension of a cell complex
function dimension(cells) {
  var d = 0
    , max = Math.max
  for(var i=0, il=cells.length; i<il; ++i) {
    d = max(d, cells[i].length)
  }
  return d-1
}
exports.dimension = dimension

//Counts the number of vertices in faces
function countVertices(cells) {
  var vc = -1
    , max = Math.max
  for(var i=0, il=cells.length; i<il; ++i) {
    var c = cells[i]
    for(var j=0, jl=c.length; j<jl; ++j) {
      vc = max(vc, c[j])
    }
  }
  return vc+1
}
exports.countVertices = countVertices

//Returns a deep copy of cells
function cloneCells(cells) {
  var ncells = new Array(cells.length)
  for(var i=0, il=cells.length; i<il; ++i) {
    ncells[i] = cells[i].slice(0)
  }
  return ncells
}
exports.cloneCells = cloneCells

//Ranks a pair of cells up to permutation
function compareCells(a, b) {
  var n = a.length
    , t = a.length - b.length
    , min = Math.min
  if(t) {
    return t
  }
  switch(n) {
    case 0:
      return 0;
    case 1:
      return a[0] - b[0];
    case 2:
      var d = a[0]+a[1]-b[0]-b[1]
      if(d) {
        return d
      }
      return min(a[0],a[1]) - min(b[0],b[1])
    case 3:
      var l1 = a[0]+a[1]
        , m1 = b[0]+b[1]
      d = l1+a[2] - (m1+b[2])
      if(d) {
        return d
      }
      var l0 = min(a[0], a[1])
        , m0 = min(b[0], b[1])
        , d  = min(l0, a[2]) - min(m0, b[2])
      if(d) {
        return d
      }
      return min(l0+a[2], l1) - min(m0+b[2], m1)
    
    //TODO: Maybe optimize n=4 as well?
    
    default:
      var as = a.slice(0)
      as.sort()
      var bs = b.slice(0)
      bs.sort()
      for(var i=0; i<n; ++i) {
        t = as[i] - bs[i]
        if(t) {
          return t
        }
      }
      return 0
  }
}
exports.compareCells = compareCells

function compareZipped(a, b) {
  return compareCells(a[0], b[0])
}

//Puts a cell complex into normal order for the purposes of findCell queries
function normalize(cells, attr) {
  if(attr) {
    var len = cells.length
    var zipped = new Array(len)
    for(var i=0; i<len; ++i) {
      zipped[i] = [cells[i], attr[i]]
    }
    zipped.sort(compareZipped)
    for(var i=0; i<len; ++i) {
      cells[i] = zipped[i][0]
      attr[i] = zipped[i][1]
    }
    return cells
  } else {
    cells.sort(compareCells)
    return cells
  }
}
exports.normalize = normalize

//Removes all duplicate cells in the complex
function unique(cells) {
  if(cells.length === 0) {
    return []
  }
  var ptr = 1
    , len = cells.length
  for(var i=1; i<len; ++i) {
    var a = cells[i]
    if(compareCells(a, cells[i-1])) {
      if(i === ptr) {
        ptr++
        continue
      }
      cells[ptr++] = a
    }
  }
  cells.length = ptr
  return cells
}
exports.unique = unique;

//Finds a cell in a normalized cell complex
function findCell(cells, c) {
  var lo = 0
    , hi = cells.length-1
    , r  = -1
  while (lo <= hi) {
    var mid = (lo + hi) >> 1
      , s   = compareCells(cells[mid], c)
    if(s <= 0) {
      if(s === 0) {
        r = mid
      }
      lo = mid + 1
    } else if(s > 0) {
      hi = mid - 1
    }
  }
  return r
}
exports.findCell = findCell;

//Builds an index for an n-cell.  This is more general than dual, but less efficient
function incidence(from_cells, to_cells) {
  var index = new Array(from_cells.length)
  for(var i=0, il=index.length; i<il; ++i) {
    index[i] = []
  }
  var b = []
  for(var i=0, n=to_cells.length; i<n; ++i) {
    var c = to_cells[i]
    var cl = c.length
    for(var k=1, kn=(1<<cl); k<kn; ++k) {
      b.length = bits.popCount(k)
      var l = 0
      for(var j=0; j<cl; ++j) {
        if(k & (1<<j)) {
          b[l++] = c[j]
        }
      }
      var idx=findCell(from_cells, b)
      if(idx < 0) {
        continue
      }
      while(true) {
        index[idx++].push(i)
        if(idx >= from_cells.length || compareCells(from_cells[idx], b) !== 0) {
          break
        }
      }
    }
  }
  return index
}
exports.incidence = incidence

//Computes the dual of the mesh.  This is basically an optimized version of buildIndex for the situation where from_cells is just the list of vertices
function dual(cells, vertex_count) {
  if(!vertex_count) {
    return incidence(unique(skeleton(cells, 0)), cells, 0)
  }
  var res = new Array(vertex_count)
  for(var i=0; i<vertex_count; ++i) {
    res[i] = []
  }
  for(var i=0, len=cells.length; i<len; ++i) {
    var c = cells[i]
    for(var j=0, cl=c.length; j<cl; ++j) {
      res[c[j]].push(i)
    }
  }
  return res
}
exports.dual = dual

//Enumerates all cells in the complex
function explode(cells) {
  var result = []
  for(var i=0, il=cells.length; i<il; ++i) {
    var c = cells[i]
      , cl = c.length|0
    for(var j=1, jl=(1<<cl); j<jl; ++j) {
      var b = []
      for(var k=0; k<cl; ++k) {
        if((j >>> k) & 1) {
          b.push(c[k])
        }
      }
      result.push(b)
    }
  }
  return normalize(result)
}
exports.explode = explode

//Enumerates all of the n-cells of a cell complex
function skeleton(cells, n) {
  if(n < 0) {
    return []
  }
  var result = []
    , k0     = (1<<(n+1))-1
  for(var i=0; i<cells.length; ++i) {
    var c = cells[i]
    for(var k=k0; k<(1<<c.length); k=bits.nextCombination(k)) {
      var b = new Array(n+1)
        , l = 0
      for(var j=0; j<c.length; ++j) {
        if(k & (1<<j)) {
          b[l++] = c[j]
        }
      }
      result.push(b)
    }
  }
  return normalize(result)
}
exports.skeleton = skeleton;

//Computes the boundary of all cells, does not remove duplicates
function boundary(cells) {
  var res = []
  for(var i=0,il=cells.length; i<il; ++i) {
    var c = cells[i]
    for(var j=0,cl=c.length; j<cl; ++j) {
      var b = new Array(c.length-1)
      for(var k=0, l=0; k<cl; ++k) {
        if(k !== j) {
          b[l++] = c[k]
        }
      }
      res.push(b)
    }
  }
  return normalize(res)
}
exports.boundary = boundary;

//Computes connected components for a dense cell complex
function connectedComponents_dense(cells, vertex_count) {
  var labels = new UnionFind(vertex_count)
  for(var i=0; i<cells.length; ++i) {
    var c = cells[i]
    for(var j=0; j<c.length; ++j) {
      for(var k=j+1; k<c.length; ++k) {
        labels.link(c[j], c[k])
      }
    }
  }
  var components = []
    , component_labels = labels.ranks
  for(var i=0; i<component_labels.length; ++i) {
    component_labels[i] = -1
  }
  for(var i=0; i<cells.length; ++i) {
    var l = labels.find(cells[i][0])
    if(component_labels[l] < 0) {
      component_labels[l] = components.length
      components.push([cells[i].slice(0)])
    } else {
      components[component_labels[l]].push(cells[i].slice(0))
    }
  }
  return components
}

//Computes connected components for a sparse graph
function connectedComponents_sparse(cells) {
  var vertices  = unique(normalize(skeleton(cells, 0)))
    , labels    = new UnionFind(vertices.length)
  for(var i=0; i<cells.length; ++i) {
    var c = cells[i]
    for(var j=0; j<c.length; ++j) {
      var vj = findCell(vertices, [c[j]])
      for(var k=j+1; k<c.length; ++k) {
        labels.link(vj, findCell(vertices, [c[k]]))
      }
    }
  }
  var components        = []
    , component_labels  = labels.ranks
  for(var i=0; i<component_labels.length; ++i) {
    component_labels[i] = -1
  }
  for(var i=0; i<cells.length; ++i) {
    var l = labels.find(findCell(vertices, [cells[i][0]]));
    if(component_labels[l] < 0) {
      component_labels[l] = components.length
      components.push([cells[i].slice(0)])
    } else {
      components[component_labels[l]].push(cells[i].slice(0))
    }
  }
  return components
}

//Computes connected components for a cell complex
function connectedComponents(cells, vertex_count) {
  if(vertex_count) {
    return connectedComponents_dense(cells, vertex_count)
  }
  return connectedComponents_sparse(cells)
}
exports.connectedComponents = connectedComponents

},{"bit-twiddle":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js","union-find":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js":[function(require,module,exports){
"use strict"

function unique_pred(list, compare) {
  var ptr = 1
    , len = list.length
    , a=list[0], b=list[0]
  for(var i=1; i<len; ++i) {
    b = a
    a = list[i]
    if(compare(a, b)) {
      if(i === ptr) {
        ptr++
        continue
      }
      list[ptr++] = a
    }
  }
  list.length = ptr
  return list
}

function unique_eq(list) {
  var ptr = 1
    , len = list.length
    , a=list[0], b = list[0]
  for(var i=1; i<len; ++i, b=a) {
    b = a
    a = list[i]
    if(a !== b) {
      if(i === ptr) {
        ptr++
        continue
      }
      list[ptr++] = a
    }
  }
  list.length = ptr
  return list
}

function unique(list, compare, sorted) {
  if(list.length === 0) {
    return list
  }
  if(compare) {
    if(!sorted) {
      list.sort(compare)
    }
    return unique_pred(list, compare)
  }
  if(!sorted) {
    list.sort()
  }
  return unique_eq(list)
}

module.exports = unique

},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/triangulate.js":[function(require,module,exports){
"use strict"

var ch = require("incremental-convex-hull")
var uniq = require("uniq")

module.exports = triangulate

function LiftedPoint(p, i) {
  this.point = p
  this.index = i
}

function compareLifted(a, b) {
  var ap = a.point
  var bp = b.point
  var d = ap.length
  for(var i=0; i<d; ++i) {
    var s = bp[i] - ap[i]
    if(s) {
      return s
    }
  }
  return 0
}

function triangulate1D(n, points, includePointAtInfinity) {
  if(n === 1) {
    if(includePointAtInfinity) {
      return [ [-1, 0] ]
    } else {
      return []
    }
  }
  var lifted = points.map(function(p, i) {
    return [ p[0], i ]
  })
  lifted.sort(function(a,b) {
    return a[0] - b[0]
  })
  var cells = new Array(n - 1)
  for(var i=1; i<n; ++i) {
    var a = lifted[i-1]
    var b = lifted[i]
    cells[i-1] = [ a[1], b[1] ]
  }
  if(includePointAtInfinity) {
    cells.push(
      [ -1, cells[0][1], ],
      [ cells[n-1][1], -1 ])
  }
  return cells
}

function triangulate(points, includePointAtInfinity) {
  var n = points.length
  if(n === 0) {
    return []
  }
  
  var d = points[0].length
  if(d < 1) {
    return []
  }

  //Special case:  For 1D we can just sort the points
  if(d === 1) {
    return triangulate1D(n, points, includePointAtInfinity)
  }
  
  //Lift points, sort
  var lifted = new Array(n)
  var upper = 1.0
  for(var i=0; i<n; ++i) {
    var p = points[i]
    var x = new Array(d+1)
    var l = 0.0
    for(var j=0; j<d; ++j) {
      var v = p[j]
      x[j] = v
      l += v * v
    }
    x[d] = l
    lifted[i] = new LiftedPoint(x, i)
    upper = Math.max(l, upper)
  }
  uniq(lifted, compareLifted)
  
  //Double points
  n = lifted.length

  //Create new list of points
  var dpoints = new Array(n + d + 1)
  var dindex = new Array(n + d + 1)

  //Add steiner points at top
  var u = (d+1) * (d+1) * upper
  var y = new Array(d+1)
  for(var i=0; i<=d; ++i) {
    y[i] = 0.0
  }
  y[d] = u

  dpoints[0] = y.slice()
  dindex[0] = -1

  for(var i=0; i<=d; ++i) {
    var x = y.slice()
    x[i] = 1
    dpoints[i+1] = x
    dindex[i+1] = -1
  }

  //Copy rest of the points over
  for(var i=0; i<n; ++i) {
    var h = lifted[i]
    dpoints[i + d + 1] = h.point
    dindex[i + d + 1] =  h.index
  }

  //Construct convex hull
  var hull = ch(dpoints, false)
  if(includePointAtInfinity) {
    hull = hull.filter(function(cell) {
      var count = 0
      for(var j=0; j<=d; ++j) {
        var v = dindex[cell[j]]
        if(v < 0) {
          if(++count >= 2) {
            return false
          }
        }
        cell[j] = v
      }
      return true
    })
  } else {
    hull = hull.filter(function(cell) {
      for(var i=0; i<=d; ++i) {
        var v = dindex[cell[i]]
        if(v < 0) {
          return false
        }
        cell[i] = v
      }
      return true
    })
  }

  if(d & 1) {
    for(var i=0; i<hull.length; ++i) {
      var h = hull[i]
      var x = h[0]
      h[0] = h[1]
      h[1] = x
    }
  }

  return hull
}
},{"incremental-convex-hull":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js","uniq":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/parse-svg-path/index.js":[function(require,module,exports){

module.exports = parse

/**
 * expected argument lengths
 * @type {Object}
 */

var length = {a: 7, c: 6, h: 1, l: 2, m: 2, q: 4, s: 4, t: 2, v: 1, z: 0}

/**
 * segment pattern
 * @type {RegExp}
 */

var segment = /([astvzqmhlc])([^astvzqmhlc]*)/ig

/**
 * parse an svg path data string. Generates an Array
 * of commands where each command is an Array of the
 * form `[command, arg1, arg2, ...]`
 *
 * @param {String} path
 * @return {Array}
 */

function parse(path) {
	var data = []
	path.replace(segment, function(_, command, args){
		var type = command.toLowerCase()
		args = parseValues(args)

		// overloaded moveTo
		if (type == 'm' && args.length > 2) {
			data.push([command].concat(args.splice(0, 2)))
			type = 'l'
			command = command == 'm' ? 'l' : 'L'
		}

		while (true) {
			if (args.length == length[type]) {
				args.unshift(command)
				return data.push(args)
			}
			if (args.length < length[type]) throw new Error('malformed path data')
			data.push([command].concat(args.splice(0, length[type])))
		}
	})
	return data
}

function parseValues(args){
	args = args.match(/-?[.0-9]+(?:e[-+]?\d+)?/ig)
	return args ? args.map(Number) : []
}

},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/ant.js":[function(require,module,exports){
'use strict';

var sign = require('./utilities.js').sign;
var calculateDistance = require('./utilities.js').distance;

// var points = require('./initializePoints.js').points;
// var citySet = require('./initializePoints.js').citySet;
// var textPointsId = require('./initializePoints.js').textPointsId;
// var possibleStartPointsId = require('./initializePoints.js').possibleStartPointsId;

var liveMousePosition = require('./mouse.js');

var Vector = require('./vector.js');

var random = Math.random;
var floor = Math.floor;
// var REPULSION = 0.05;
// var REPULSIONSPEED = 0.002;
// var ANTVELOCITY = 0.001;

module.exports = function(container, initPoints, options){

    console.log('Options ant :', options);
    // Define those parameters as attributes of Ant object ?
    var REPULSION = options.repSize;
    var REPULSIONSPEED = options.repSpeed;
    var ANTVELOCITY = options.velocity;
    var WEIGHT = options.weight;

    var mouse = liveMousePosition(container);

    var points = initPoints.points;
    var citySet = initPoints.citySet;
    var textPointsId = initPoints.textPointsId;
    var possibleStartPointsId = initPoints.possibleStartPointsId;


    function Ant(point) {
        this.x = point.x;                
        this.y = point.y;
        this.velocity = ANTVELOCITY;
        this.weight = WEIGHT;
        this.repSize = REPULSION;
        this.repSpeed = REPULSIONSPEED;
        this.edge = undefined;
        this.state = "forage";
        this.edges = [];
        this.lastCity = undefined;
        this.origin = point;
        this.destination = undefined;
        this.orientation = undefined;
        this.direction = new Vector(0,0);
        this.prog = 0;
    }
    // forage: the ant wanders around without any pheromon deposition
    // once it finds a city, it starts remembering the nodes it goes through
    // when it finds another city, it computes the path length and adds pheromons one each edges
    // proportionnaly to the shortestness of the path
    // it resets the list of nodes and continues
    // while foraging the ant choses the path with a pheromon preference


    // static methods
    Ant.generateRandStartPoint = function() {
        var randId = Math.floor(possibleStartPointsId.length * random());
        var randStartPoint = points[possibleStartPointsId[randId]];
        return randStartPoint;
    }


    // methods
    Ant.prototype = {

        transit: function(){
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
                        var antWeight = this.weight;
                        this.edges.forEach(function(e){
                            var a = e.pt1, b = e.pt2, weight = 1;  
                            // increased dropped pheromons for textEdges
                            if ((citySet.indexOf(a.id) != -1) && citySet.indexOf(b.id) != -1 && (Math.abs(a.id - b.id) == 1))
                            {
                                weight *= antWeight;
                            }
                            e.pheromon += (deltaPheromone * weight);
                        });

                        this.edges = [this.edge];
                        this.lastCity = this.origin.id;
                    }
                }
              break;
            }

        },

        setDirection: function(){
            var possibleEdges = [];

            for (var i = 0; i < this.origin.nexts.length; i++)
            {
                possibleEdges[i] = this.origin.nexts[i];
            } 

            possibleEdges.splice(possibleEdges.indexOf(this.edge),1);

            // flip a coin and either take the smelliest path or a random one
            if (random() > 0.5){
                var smells = possibleEdges.map(function(e){return e.pheromon});
                var index = smells.indexOf(Math.max.apply(Math, smells));
                this.edge = possibleEdges[index];
            } 
            else
                this.edge = possibleEdges[floor(random()*possibleEdges.length)];

            // set the destination point, being edge.pt1 or edge.pt2
            this.destination = (this.origin == this.edge.pt1) ? this.edge.pt2 : this.edge.pt1;

            this.direction.x = this.destination.x - this.origin.x; 
            this.direction.y = this.destination.y - this.origin.y;

            this.direction.normalize();
        },

        move: function(){
            var edgeChanged;
            var cityReached = false;

            this.direction.x = this.destination.x - this.x; 
            this.direction.y = this.destination.y - this.y;
            this.direction.normalize();

            // on edge
            if ((calculateDistance(this, this.destination) > this.repSpeed)){

                // a delta movement will be applied if collision with obstacle detected
                var delta = this.avoidObstacle();

                this.x += this.velocity * this.direction.x + delta.x * this.repSpeed;
                this.y += this.velocity * this.direction.y + delta.y * this.repSpeed;

                this.prog = this.calculateProgression();
                
                edgeChanged = false;

            // on vertex
            } else {
                this.step = 0;
                this.prog = 0;
                this.origin = this.destination;
                this.x = this.origin.x;
                this.y = this.origin.y;

                this.setDirection();

                cityReached = (citySet.indexOf(this.origin.id) != -1);
                edgeChanged = true;
            }
            return {cityReached: cityReached, edgeChanged: edgeChanged};
        },

        avoidObstacle: function(){
            var distance = calculateDistance(this, mouse);
        
            if (distance <= this.repSize) {

                return {
                    // delta movement is composed of a repulsion delta and a circular delta 
                    x: (this.x - mouse.x)/distance + (this.y - mouse.y)/distance * 1,
                    y: (this.y - mouse.y)/distance - (this.x - mouse.x)/distance * 1
                };
            }
            else
                return {x:0, y:0};
        },

        calculateProgression: function(){
            var v = new Vector(this.x - this.origin.x, this.y - this.origin.y);
            var norm = v.norm();

            var theta = (v.x * this.edge.direction.x + v.y * this.edge.direction.y) / norm;
            var prog = norm * Math.abs(theta);
            // returns length of projection on edge
            return prog;
        }

    };
    return Ant;
}


},{"./mouse.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/mouse.js","./utilities.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/utilities.js","./vector.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/vector.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/createEdges.js":[function(require,module,exports){
'use strict'

var dt = require("delaunay-triangulate");

var range = require('./utilities.js').range;

// var points = require('./initializePoints.js').points;
var textMesh = require('./initializePoints.js').textMesh;
var citySet = require('./initializePoints.js').citySet;
var nbRandomPoints = require('./initializePoints.js').nbRandomPoints;
var forcedEdges = require('./initializePoints.js').forcedEdges;

var Edge = require('./edge.js');


module.exports = function(points){
    // triangulate
    var cells = dt(points.map(function(p){
        return [p.x, p.y];
    }));

    var edges = [];
    var permutations = [[0,1], [0,2], [1,2]];

    // force the edges of the text to be edges of the graph
    if (textMesh) {
        range(0, points.length - nbRandomPoints).forEach(function(id){
            var directLink = forcedEdges[id];
            var textEdge = Edge.create(points[id], points[directLink]);
            edges.push(textEdge);
            points[id].nexts.push(textEdge);
        })
    }


    cells.forEach(function(cell){
       
        for (var i = 0; i < 3; ++i){  // for each point.id listed in current cell
            var pt = points[cell[i]];

            for (var j = 1; j < 3; ++j){ 

                var ptj = points[cell[( i + j ) % 3]]; // pick one of the other 2 points of the cell
                var newEdge = undefined;

                // if pt already has nextEdges ...
                if (pt.nexts.length != 0) {
                    
                    // ... get the points corresponding ...
                    var tempPoints = pt.nexts.map(function(e){
                        return [e.pt1, e.pt2];
                    }).reduce(function(a, b){
                         return a.concat(b);
                    });

                    // ... and check if ptj already is part of the existing nextEdges. If not, add the edge.
                    if (tempPoints.indexOf(ptj) == -1){
                        newEdge = Edge.create(pt, ptj);
                        edges.push(newEdge);
                        pt.nexts.push(newEdge);
                    }
                }
                else {
                    newEdge = Edge.create(pt, ptj);
                    edges.push(newEdge);
                    pt.nexts.push(newEdge);
                }

                // add also the edge to the edge's other point's nextEdges
                if (newEdge != undefined){
                    ptj.nexts.push(newEdge);
                }         
            }

            // add the textEdges to nextEdges map
            if (textMesh && (citySet.indexOf(pt) != -1)) {
                var textEdge = Edge.create(pt, points[pt.id + 1]);
                edges.push(textEdge);
                pt.nexts.push(textEdge);
            }

        }
    });

    return edges;
};
},{"./edge.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/edge.js","./initializePoints.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/initializePoints.js","./utilities.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/utilities.js","delaunay-triangulate":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/delaunay-triangulate/triangulate.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/edge.js":[function(require,module,exports){
'use strict';

var sqrt = Math.sqrt;
var pow = Math.pow;
var abs = Math.abs;
var atan = Math.atan;

var Vector = require('./vector.js');


function Edge(ptA, ptB) {
    var distance = sqrt( pow(ptA.x - ptB.x, 2) + pow(ptA.y - ptB.y, 2) );

    // find line equation ax + by + c = 0
    var a = 1;
    var b = - (ptB.x - ptA.x) / (ptB.y - ptA.y);

    // orientate vector (a,b)
    if (b < 0){
        b = -b;
        a = -a;
    }

    // normalize vector (a,b)
    var n = new Vector(a, b);
    n.normalize();

    var c = - (a * ptA.x + b * ptA.y);

    // // calculate vector director
    var v = new Vector(ptB.x - ptA.x, ptB.y - ptA.y);
    
    v.normalize();


    this.id = undefined;
    this.pt1 = ptA;
    this.pt2 = ptB;
    this.direction = v;
    this.orthDirection = n; 
    this.distance = distance;
    this.pheromon = 1/distance;
    this.line = {
        a: a,
        b: b,
        c: c,
    }
}


// static methods
Edge.create = function(ptA, ptB) {
    var edge = new Edge(ptA, ptB);
    return edge;
}


// methods
Edge.prototype = {

    getOtherPoint: function(point) {
        if (point == this.pt1)
            return this.pt2;
        else if (point == this.pt2)
            return this.pt1;
        else
            console.log("Error");
    },

    calculateDistance: function(x, y) {
        var a = this.line.a,
            b = this.line.b,
            c = this.line.c;
        return abs(a * x + b * y + c) / Math.sqrt(Math.pow(a,2) + Math.pow(b,2));
    },

}
module.exports = Edge;
},{"./vector.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/vector.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/initializeAnts.js":[function(require,module,exports){
'use strict'

var antFunction = require('./ant.js');

// var NBANTS = 4000;

module.exports = function (container, pointsInfos, options) {

	var Ant = antFunction(container, pointsInfos, options);

	var population = new Array(options.nbAnts);
	var possibleStartPointsId = pointsInfos.possibleStartPointsId;

	for (var i = 0; i < options.nbAnts; i++) {
		var newAnt = new Ant(Ant.generateRandStartPoint());
		newAnt.setDirection();
		population[i] = newAnt;
	}

	return population;

}
},{"./ant.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/ant.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/initializePoints.js":[function(require,module,exports){
'use strict'

var parse = require('parse-svg-path');

var range = require('./utilities.js').range;

var Point = require('./point.js');

var random = Math.random;

// var nbRandomPoints = 500;
// var nbStartPoints = 20;

var nbCity = 2;

var textMesh = true;

// Frame definition
var xInit = 0, yInit = 0;
var w = 1,
    h = 1;

var Achar = "c 4.6011,-11.71047 9.20835,-23.42006 13.8199,-35.12898 4.61156,-11.70892 9.22741,-23.41718 13.84573,-35.125 4.61831,-11.70782 9.23908,-23.41519 13.86046,-35.12233 4.62138,-11.70714 9.24336,-23.41406 13.86411,-35.12097 4.62074,-11.70691 9.24025,-23.41381 13.85667,-35.12092 4.61641,-11.70712 9.22974,-23.41444 13.83814,-35.12218 4.60839,-11.70775 9.21186,-23.41592 13.80854,-35.12474 4.59668,-11.70881 9.18658,-23.418277 13.76785,-35.128603 12.99923,-3.357855 24.30069,-6.821144 33.86893,-4.46411 9.56825,2.357035 17.40328,10.534393 23.46967,30.457833 4.17598,10.62114 8.36031,21.23882 12.54942,31.85453 4.1891,10.61571 8.38298,21.22943 12.57807,31.84266 4.19509,10.61323 8.3914,21.22595 12.58535,31.83965 4.19395,10.61369 8.38555,21.22836 12.57123,31.84549 4.18568,10.61712 8.36544,21.2367 12.53572,31.86021 4.17028,10.6235 8.33108,21.25093 12.47883,31.88377 4.14775,10.63283 8.28245,21.27107 12.40053,31.9162 4.11808,10.64512 8.21955,21.29712 12.30085,31.95749 -5.90896,6.95561 -24.61617,1.11298 -35.59372,3 -16.75248,2.84859 -26.96421,-0.41416 -28.40628,-19 -3.17726,-8.42157 -6.35606,-16.87924 -9.47997,-25.34854 -3.12391,-8.46929 -6.19293,-16.9502 -9.15062,-25.41826 -16.46723,-0.80053 -35.17768,-2.74281 -53.33727,-3.11439 -18.1596,-0.37158 -35.76834,0.82753 -50.03214,6.30976 -4.15679,10.93124 -8.13806,21.9269 -12.15785,32.90546 -4.0198,10.97855 -8.07813,21.94001 -12.38903,32.80282 -9.52758,0.82747 -19.10457,0.96423 -28.69282,0.93363 -9.58824,-0.0306 -19.18773,-0.22855 -28.7603,-0.0705 l 0,-2.19791 z m 174,-109 c -3.17462,-9.51136 -6.44835,-18.99174 -9.7543,-28.46224 -6.81559,-19.51412 -6.24202,-25.37946 -12.3174,-43.60788 -3.16722,-9.51543 -13.60257,-32.32866 -16.49973,-41.92988 -6.08989,12.88576 -11.132,26.59272 -15.93415,40.47476 -4.80215,13.88203 -9.36435,27.93915 -14.49442,41.52524 -2.39425,12.05801 -22.8815,40.11116 2.18745,34 11.05256,-0.44486 22.52939,0.11061 33.85624,0.24955 11.32684,0.13895 22.50369,-0.13863 32.95631,-2.24955 z ";
var Nchar = "c 0,0 0,-23.83334 0,-35.75 0,-11.91667 0,-23.83334 0,-35.75 0,-11.91666 0,-23.83333 0,-35.75 0,-11.916667 0,-23.833335 0,-35.750003 14.64729,1.05313 31.03193,-2.08808 44.60679,1.55415 6.18076,7.610996 12.35438,15.231342 18.51973,22.861013 6.16535,7.62967 12.32244,15.26866 18.47014,22.91696 6.1477,7.64829 12.28602,15.30588 18.41383,22.97274 6.12781,7.66686 12.24512,15.343 18.35081,23.02838 6.10568,7.68538 12.19975,15.38001 18.28107,23.08386 6.08132,7.70384 12.1499,15.41691 18.20462,23.13918 6.05472,7.72226 12.09557,15.45372 18.12145,23.19435 6.02587,7.74063 12.03676,15.49043 18.03156,23.24937 0.4558,-15.4914 0.72655,-30.98725 0.88119,-46.48582 0.15464,-15.49858 0.19319,-30.99988 0.18458,-46.50218 -0.009,-15.5023 -0.0643,-31.0056 -0.0983,-46.50817 -0.034,-15.50258 -0.0461,-31.004432 0.0325,-46.503833 9.33333,0 18.66667,0 28,0 9.33333,0 18.66666,0 28,0 0,11.916667 0,23.833332 0,35.750003 0,11.91666 0,23.83333 0,35.75 0,11.91666 0,23.83333 0,35.75 0,11.91667 0,23.83333 0,35.75 0,11.91666 0,23.83333 0,35.75 0,11.91667 0,23.83334 0,35.75 0,11.91666 0,23.83333 0,35.75 0,11.91667 0,23.83334 0,35.75 -14.6441,-1.05348 -31.026,2.08847 -44.59764,-1.55416 -6.27247,-7.38037 -12.48436,-14.81692 -18.65123,-22.29507 -6.16688,-7.47815 -12.28873,-14.9979 -18.38111,-22.54465 -6.09238,-7.54675 -12.15528,-15.12051 -18.20425,-22.70667 -6.04896,-7.58617 -12.08399,-15.18476 -18.12063,-22.78116 -6.03664,-7.5964 -12.07489,-15.19062 -18.1303,-22.76807 -6.05541,-7.57745 -12.12797,-15.13813 -18.23323,-22.66744 -6.10526,-7.52932 -12.24322,-15.02727 -18.42942,-22.47926 -6.18619,-7.45199 -12.42063,-14.85803 -18.71886,-22.20352 -0.22791,15.16496 -0.36674,30.3308 -0.4489,45.49718 -0.0822,15.16637 -0.10766,30.33329 -0.10889,45.5004 -0.001,15.16712 0.0218,30.33443 0.0367,45.50162 0.0149,15.16718 0.0216,30.33422 -0.0122,45.5008 -9.33333,0 -18.66667,0 -28,0 -9.33333,0 -18.66666,0 -28,0 0,-11.91667 0,-23.83334 0,-35.75 0,-11.91667 0,-23.83334 0,-35.75 0,-11.91666 0,-23.83333 0,-35.75 0,-11.91667 0,-35.75 -10e-6,-35.75 z ";
var Tchar = "c 0,-9.91667 0,-19.83334 0,-29.75 0,-9.91667 0,-19.83334 0,-29.75 0,-9.91666 0,-19.83333 0,-29.75 0,-9.91666 -1.47394,-19.83333 -1.47394,-29.75 0,-15.33334 -29.19273,0 -44.52606,0 -15.33333,0 -30.66667,0 -46,0 0,-8 0,-16 0,-24.000002 0,-8 0,-16.000001 0,-24.000001 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 9.91666,0 19.83333,0 29.75,0 9.91667,0 19.83333,0 29.75,0 0,8 0,16.000001 0,24.000001 0,8.000002 0,16.000002 0,24.000002 -15,0 -30,0 -45,0 -15,0 -44.88809,-13.52565 -45,1.47394 -0.0735,9.85588 -0.10613,18.2399 -0.11249,28.09923 -0.006,9.85932 0.0135,19.72001 0.045,29.58133 0.0315,9.86132 0.0745,19.72329 0.1143,29.58517 0.0399,9.86188 0.0765,19.72367 0.0954,29.58467 0.0189,9.86099 0.0199,19.72118 -0.0117,29.57984 -0.0315,9.85866 -0.0956,19.7158 -0.20688,29.57069 -0.11131,9.85488 -0.26985,19.70752 -0.49033,29.55719 -0.22047,9.84967 -0.50289,19.69636 -0.86193,29.53937 -7.39624,-0.99964 -18.95658,0.96726 -29.70912,1.82971 -10.75253,0.86245 -24.48556,2.02609 -24.86231,-4.79696 -0.66223,-11.99314 -0.66223,-21.54349 -0.49667,-30.48314 0.16556,-8.93965 0.49667,-17.2686 0.49667,-26.81895 0,-9.55035 0,-19.1007 0,-28.65105 0,-9.55035 0,-19.10069 0,-28.65104 z";
var Schar = "c -16.6587,-2.34387 -33.35995,-6.23515 -49.25137,-12.08831 -15.89142,-5.85316 -30.97305,-13.6682 -44.3926,-23.85957 3.73715,-8.02639 7.96065,-15.79371 12.19101,-23.55228 4.23037,-7.75858 8.46759,-15.50842 12.23219,-23.49984 8.82902,6.69162 18.48235,12.75663 28.67412,17.86344 10.19178,5.10681 20.92205,9.25542 31.90485,12.11421 10.9829,2.8588 22.2183,4.4278 33.4206,4.37539 11.2022,-0.0524 22.3713,-1.72622 33.2212,-5.35304 14.2804,-6.92169 17.3033,-19.64436 13.4946,-31.15028 -3.8087,-11.50591 -14.4489,-21.79507 -27.4946,-23.84972 -10.162,-4.07686 -21.1052,-7.14863 -32.1975,-10.11097 -11.0924,-2.96235 -22.334,-5.81526 -33.0931,-9.45439 -10.7591,-3.63914 -21.03562,-8.06449 -30.19779,-14.17171 -9.16216,-6.10723 -17.20996,-13.89632 -23.51158,-24.26293 -5.1794,-10.97058 -7.48544,-22.77583 -7.30512,-34.52966 0.18031,-11.75382 2.84698,-23.4562 7.61299,-34.22103 4.76601,-10.76483 11.63137,-20.5921 20.20905,-28.595692 8.57769,-8.003592 18.86771,-14.183506 30.48305,-17.653621 12.1821,-4.024145 24.8777,-6.357009 37.6988,-7.121901 12.8211,-0.764892 25.7677,0.03819 38.4517,2.285933 12.6841,2.247744 25.1055,5.940153 36.8765,10.953918 11.7709,5.013764 22.8912,11.348885 32.973,18.882053 -4.4745,7.00368 -8.2581,14.70708 -12.1886,22.23499 -3.9306,7.5279 -8.008,14.88031 -13.0701,21.18198 -15.0728,-11.15158 -33.5009,-20.54491 -52.6114,-24.91726 -19.1105,-4.37235 -38.9034,-3.72372 -56.7058,5.20861 -10.0321,7.02789 -13.3348,18.84695 -11.1561,29.53596 2.1787,10.68902 9.8387,20.24799 21.732,22.75572 10.1335,4.43307 21.0241,7.69407 32.1101,10.69023 11.0861,2.99616 22.3676,5.72747 33.283,9.10117 10.9155,3.3737 21.4648,7.38978 31.0865,12.95547 9.6217,5.56569 18.3157,12.68099 25.5204,22.25313 6.6301,9.71635 10.5587,20.88023 12.0064,32.41231 1.4476,11.53207 0.4142,23.43234 -2.8797,34.62146 -3.2939,11.18912 -8.8484,21.66709 -16.443,30.35457 -7.5946,8.68749 -17.2293,15.58449 -28.6837,19.61166 -13.0091,5.92928 -27.0197,8.63869 -41.2728,9.63608 -14.253,0.99738 -28.7484,0.28274 -42.7272,-0.63608 z ";

var svgString = "m 0,0 " + Achar +"m 126,-31 " + Nchar + "m 376,24 "+ Tchar +"m 250,120 "+ Schar;

 
function svgToPoints(svgString) {
    var points = [];
    var edges = Object.create(null);

    var beginingPath;

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
                beginingPath = nbPoints;
                break;
            case "M":
                X = command[1];
                Y = command[2];
                prevPoint = undefined;
                beginingPath = nbPoints;
                break;  
            case "c":
                X += command[5];
                Y += command[6];
                points.push({id:nbPoints, x:X, y:Y});

                if (prevPoint != undefined) {
                    edges[prevPoint] = nbPoints;
                }
                prevPoint = nbPoints;
                nbPoints++;
                break;
            case "z":
                edges[prevPoint] = nbPoints;
                beginingPath = undefined;
                prevPoint = undefined;
                break;    
        }
    }
    return {points : points, edges : edges};
}

// initialize points

module.exports = function(nbStartPoints, nbRandomPoints){
    var points = [];
    var forcedEdges;
    var citySet;

    if (textMesh){

        var myText = svgToPoints(svgString);
        points = myText.points;
        forcedEdges = myText.edges;
        citySet = range(0, points.length);

        var scaleX = 0.5;
        var scaleY = 0.2;
        var deltaX = 0.25;
        var deltaY = 0.3;

        // scale points to [0,1] + delta
        var maxX = Math.max.apply(Math, points.map(function(p){return p.x}));
        var minX = Math.min.apply(Math, points.map(function(p){return p.x}));
        var maxY = Math.max.apply(Math, points.map(function(p){return p.y}));
        var minY = Math.min.apply(Math, points.map(function(p){return p.y}));
        points = points.map(function(p){
            var x = scaleX * (p.x-minX)/(maxX-minX) + deltaX;
            var y = scaleY * (p.y-minY)/(maxY-minY) + deltaY;
            var newPoint = new Point(x, y);
            newPoint.id = p.id;

            return newPoint;
        });

        // only add random points
        var nbPoints = points.length;
        for(var i=0; i<nbRandomPoints; ++i) {

            var x = random();
            var y = random();

            var newPoint = new Point(x, y);
            newPoint.id = nbPoints;

            points.push(newPoint);

            nbPoints++;
        }

    } else {
        //add random points

        var nbPoints = 0;
        for(var i=0; i<nbRandomPoints; ++i) {

            var x = random();
            var y = random();

            var newPoint = new Point(x, y);
            newPoint.id = nbPoints;

            points.push(newPoint);
            
            nbPoints++;
        }

        citySet = range(0, nbCity);
        console.log(citySet);
    }


    // initialize start points
    var possibleStartPointsId = [];

    if (nbRandomPoints === 0) {
        for (var i = 0; i < nbStartPoints; i++){
            possibleStartPointsId.push(Math.floor(nbRandomPoints * random()));
        }
    } else {
        for (var i = 0; i < nbStartPoints; i++){
            possibleStartPointsId.push(Math.floor(nbRandomPoints * random()));
        }
    }
    

    return {
        textMesh: textMesh,
        points: points,
        citySet: citySet,
        possibleStartPointsId: possibleStartPointsId,
        nbRandomPoints: nbRandomPoints,
        forcedEdges: forcedEdges
    };
}

},{"./point.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/point.js","./utilities.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/utilities.js","parse-svg-path":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/node_modules/parse-svg-path/index.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/mouse.js":[function(require,module,exports){
'use strict'

module.exports = function (container){

	var mouse = {
	    x: 0,
	    y: 0
	};

	container.addEventListener( 'mousemove', function(e){
	    var rect = container.getBoundingClientRect();

	    mouse.x = (e.clientX - rect.left ) / rect.width;
	    mouse.y = (e.clientY - rect.top )/ rect.height;
	});

	return mouse;

};

},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/point.js":[function(require,module,exports){
'use strict'

function Point(x, y) {
    this.id = undefined;                
    this.x = x;
    this.y = y;
    this.nexts = [];
}

module.exports = Point;
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/rendering.js":[function(require,module,exports){
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

	var animID;

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

	function reset(){
		population = [];
		edges = [];
		pointsInfos = [];

		cancelAnimationFrame(animID);
	}
	
	// container.addEventListener('click', togglePlayPause);

	function animate(){
		tick();
		
		if(!paused)
			animID = requestAnimationFrame(animate);
	}
	animate();


	function setAntCount(opts){

		var previousCount = population.length;

		console.log('nbAnts :', opts.nbAnts);
		console.log('previousCount :', previousCount);


		if (opts.nbAnts > previousCount){
			console.log('pointsInfos :', pointsInfos);
			var Ant = antFunction(container, pointsInfos, opts);

			for (var i = 0; i < opts.nbAnts - previousCount; i++) {
				console.log('test');
				var newAnt = new Ant(Ant.generateRandStartPoint());
				newAnt.setDirection();
				population.push(newAnt);
			}
		}
		else{
			population = population.slice(0, opts.nbAnts);
			nbAnts = population.length;
			console.log('Nb Ants :', population.length);
		}
		console.log('Nouveau compte :', population.length);
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
		reset: reset,
		// should be a getter/setter, but IE8
		getAntCount: function(){
			return population.length;
		},
		setAntCount: setAntCount,
		modifyAnts: modifyAnts
	}
}

},{"./ant.js":"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/ant.js"}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/utilities.js":[function(require,module,exports){
'use strict';

var sqrt = Math.sqrt;
var pow = Math.pow;

function sign(x) {
	return x ? x < 0 ? -1 : 1 : 0;
}

function range(start, count) {
    return Array.apply(0, Array(count)).map(function (element, index) {
    	return index + start
    });
}

function distance(a, b){
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

function norm(v){
	return sqrt(pow(v.x, 2) + pow(v.y, 2));
}

module.exports = {
	sign: sign,
	range: range,
	distance: distance,
	norm: norm
}
},{}],"/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/src/vector.js":[function(require,module,exports){
'use strict'

function Vector(x, y) {
    this.x = x;                
    this.y = y;
}

Vector.prototype.norm = function(){
	return Math.sqrt(this.x * this.x + this.y * this.y);
}

Vector.prototype.normalize = function(){
	var norm = this.norm();
	this.x = this.x / norm;
	this.y = this.y / norm;
}



module.exports = Vector;
},{}]},{},["/Users/Romain/Documents/Programmation/Projets_Ants/AntColony/index.js"])
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi91c3IvbG9jYWwvbGliL25vZGVfbW9kdWxlcy93YXRjaGlmeS9ub2RlX21vZHVsZXMvYnJvd3NlcmlmeS9ub2RlX21vZHVsZXMvYnJvd3Nlci1wYWNrL19wcmVsdWRlLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L2luZGV4LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvaWNoLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL25vZGVfbW9kdWxlcy90d28tc3VtL3R3by1zdW0uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL25vZGVfbW9kdWxlcy9yb2J1c3Qtc2NhbGUvcm9idXN0LXNjYWxlLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1YnRyYWN0L3JvYnVzdC1kaWZmLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1bS9yb2J1c3Qtc3VtLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvdHdvLXByb2R1Y3QvdHdvLXByb2R1Y3QuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL29yaWVudGF0aW9uLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3NpbXBsaWNpYWwtY29tcGxleC9ub2RlX21vZHVsZXMvYml0LXR3aWRkbGUvdHdpZGRsZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19BbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvbm9kZV9tb2R1bGVzL3VuaW9uLWZpbmQvaW5kZXguanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvc2ltcGxpY2lhbC1jb21wbGV4L3RvcG9sb2d5LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvdW5pcS91bmlxLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS90cmlhbmd1bGF0ZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19BbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvcGFyc2Utc3ZnLXBhdGgvaW5kZXguanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvc3JjL2FudC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19BbnRzL0FudENvbG9ueS9zcmMvY3JlYXRlRWRnZXMuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvc3JjL2VkZ2UuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvc3JjL2luaXRpYWxpemVBbnRzLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L3NyYy9pbml0aWFsaXplUG9pbnRzLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L3NyYy9tb3VzZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19BbnRzL0FudENvbG9ueS9zcmMvcG9pbnQuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL1Byb2pldHNfQW50cy9BbnRDb2xvbnkvc3JjL3JlbmRlcmluZy5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vUHJvamV0c19BbnRzL0FudENvbG9ueS9zcmMvdXRpbGl0aWVzLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9Qcm9qZXRzX0FudHMvQW50Q29sb255L3NyYy92ZWN0b3IuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2pEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzdiQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2pEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDM0pBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMzSkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM3TEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNU1BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3pEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN0VkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDekRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM5SkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2REE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM3TUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNyRkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzdFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNyQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDM0tBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDbkJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ1RBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3pNQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzVCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gZSh0LG4scil7ZnVuY3Rpb24gcyhvLHUpe2lmKCFuW29dKXtpZighdFtvXSl7dmFyIGE9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtpZighdSYmYSlyZXR1cm4gYShvLCEwKTtpZihpKXJldHVybiBpKG8sITApO3ZhciBmPW5ldyBFcnJvcihcIkNhbm5vdCBmaW5kIG1vZHVsZSAnXCIrbytcIidcIik7dGhyb3cgZi5jb2RlPVwiTU9EVUxFX05PVF9GT1VORFwiLGZ9dmFyIGw9bltvXT17ZXhwb3J0czp7fX07dFtvXVswXS5jYWxsKGwuZXhwb3J0cyxmdW5jdGlvbihlKXt2YXIgbj10W29dWzFdW2VdO3JldHVybiBzKG4/bjplKX0sbCxsLmV4cG9ydHMsZSx0LG4scil9cmV0dXJuIG5bb10uZXhwb3J0c312YXIgaT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2Zvcih2YXIgbz0wO288ci5sZW5ndGg7bysrKXMocltvXSk7cmV0dXJuIHN9KSIsIid1c2Ugc3RyaWN0JztcblxudmFyIGluaXRSZW5kZXJpbmcgPSByZXF1aXJlKCcuL3NyYy9yZW5kZXJpbmcuanMnKTtcbnZhciBpbml0aWFsaXplUG9pbnRzID0gcmVxdWlyZSgnLi9zcmMvaW5pdGlhbGl6ZVBvaW50cy5qcycpO1xudmFyIGNyZWF0ZUVkZ2VzID0gcmVxdWlyZSgnLi9zcmMvY3JlYXRlRWRnZXMuanMnKTtcbnZhciBpbml0QW50cyA9IHJlcXVpcmUoJy4vc3JjL2luaXRpYWxpemVBbnRzJyk7XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24gaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRpb25zKXtcblxuXHR2YXIgcmVuZGVyLCBwb2ludHNJbmZvcywgZWRnZXMsIHBvcHVsYXRpb24sIGluaXRWYXI7XG5cblxuXHRmdW5jdGlvbiBfaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRpb25zKXtcblx0XHRwb2ludHNJbmZvcyA9IGluaXRpYWxpemVQb2ludHMob3B0aW9ucy5uYlN0YXJ0LCBvcHRpb25zLm5iUmFuZCk7XG5cdFx0ZWRnZXMgPSBjcmVhdGVFZGdlcyhwb2ludHNJbmZvcy5wb2ludHMpO1xuXHRcdHBvcHVsYXRpb24gPSBpbml0QW50cyhjb250YWluZXJFbGVtZW50LCBwb2ludHNJbmZvcywgb3B0aW9ucyk7XG5cdFx0aW5pdFZhciA9IHtcblx0XHRcdHBvaW50c0luZm9zOiBwb2ludHNJbmZvcyxcblx0XHRcdGVkZ2VzOiBlZGdlcyxcblx0XHRcdHBvcHVsYXRpb246IHBvcHVsYXRpb25cblx0XHR9O1xuXHRcdHJlbmRlciA9IGluaXRSZW5kZXJpbmcoY29udGFpbmVyRWxlbWVudCwgaW5pdFZhcik7XG5cdH1cblxuXHRfaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRpb25zKTtcblxuXHRyZXR1cm4ge1xuXHRcdHRvZ2dsZVBsYXlQYXVzZTogZnVuY3Rpb24oKXsgcmVuZGVyLnRvZ2dsZVBsYXlQYXVzZSgpIH0sXG5cdFx0Y2hhbmdlT3B0aW9uczogZnVuY3Rpb24ob3B0cyl7XG5cdFx0XHRpZiAob3B0cy5uYlN0YXJ0ID09PSB1bmRlZmluZWQgJiYgb3B0cy5uYlJhbmQgPT09IHVuZGVmaW5lZCl7XG5cdFx0XHRcdC8vIG1vZGlmeSBwb3B1bGF0aW9uXG5cdFx0XHRcdHJlbmRlci5tb2RpZnlBbnRzKG9wdHMpO1xuXHRcdFx0fVxuXHRcdFx0ZWxzZXtcblx0XHRcdFx0Ly8gcmVzZXQgcHJldmlvdXMgYW5pbWF0aW9uXG5cdFx0XHRcdHJlbmRlci5yZXNldCgpO1xuXG5cdFx0XHRcdC8vIHJlc2V0IGVsZW1lbnRzXG5cdFx0XHRcdF9pbml0KGNvbnRhaW5lckVsZW1lbnQsIG9wdHMpO1xuXHRcdFx0XHRcblx0XHRcdH1cblx0XHR9LFxuXHRcdHJlc2V0OiBmdW5jdGlvbihvcHRzKXtcblx0XHRcdHJlbmRlci5yZXNldCgpO1xuXG5cdFx0XHRcdC8vIHJlc2V0IGVsZW1lbnRzXG5cdFx0XHRfaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRzKTtcblx0XHR9XG5cdH07XG59OyIsIlwidXNlIHN0cmljdFwiXG5cbi8vSGlnaCBsZXZlbCBpZGVhOlxuLy8gMS4gVXNlIENsYXJrc29uJ3MgaW5jcmVtZW50YWwgY29uc3RydWN0aW9uIHRvIGZpbmQgY29udmV4IGh1bGxcbi8vIDIuIFBvaW50IGxvY2F0aW9uIGluIHRyaWFuZ3VsYXRpb24gYnkganVtcCBhbmQgd2Fsa1xuXG5tb2R1bGUuZXhwb3J0cyA9IGluY3JlbWVudGFsQ29udmV4SHVsbFxuXG52YXIgb3JpZW50ID0gcmVxdWlyZShcInJvYnVzdC1vcmllbnRhdGlvblwiKVxudmFyIGNvbXBhcmVDZWxsID0gcmVxdWlyZShcInNpbXBsaWNpYWwtY29tcGxleFwiKS5jb21wYXJlQ2VsbHNcblxuZnVuY3Rpb24gY29tcGFyZUludChhLCBiKSB7XG4gIHJldHVybiBhIC0gYlxufVxuXG5mdW5jdGlvbiBTaW1wbGV4KHZlcnRpY2VzLCBhZGphY2VudCwgYm91bmRhcnkpIHtcbiAgdGhpcy52ZXJ0aWNlcyA9IHZlcnRpY2VzXG4gIHRoaXMuYWRqYWNlbnQgPSBhZGphY2VudFxuICB0aGlzLmJvdW5kYXJ5ID0gYm91bmRhcnlcbiAgdGhpcy5sYXN0VmlzaXRlZCA9IC0xXG59XG5cblNpbXBsZXgucHJvdG90eXBlLmZsaXAgPSBmdW5jdGlvbigpIHtcbiAgdmFyIHQgPSB0aGlzLnZlcnRpY2VzWzBdXG4gIHRoaXMudmVydGljZXNbMF0gPSB0aGlzLnZlcnRpY2VzWzFdXG4gIHRoaXMudmVydGljZXNbMV0gPSB0XG4gIHZhciB1ID0gdGhpcy5hZGphY2VudFswXVxuICB0aGlzLmFkamFjZW50WzBdID0gdGhpcy5hZGphY2VudFsxXVxuICB0aGlzLmFkamFjZW50WzFdID0gdVxufVxuXG5mdW5jdGlvbiBHbHVlRmFjZXQodmVydGljZXMsIGNlbGwsIGluZGV4KSB7XG4gIHRoaXMudmVydGljZXMgPSB2ZXJ0aWNlc1xuICB0aGlzLmNlbGwgPSBjZWxsXG4gIHRoaXMuaW5kZXggPSBpbmRleFxufVxuXG5mdW5jdGlvbiBjb21wYXJlR2x1ZShhLCBiKSB7XG4gIHJldHVybiBjb21wYXJlQ2VsbChhLnZlcnRpY2VzLCBiLnZlcnRpY2VzKVxufVxuXG5mdW5jdGlvbiBiYWtlT3JpZW50KGQpIHtcbiAgdmFyIGNvZGUgPSBbXCJmdW5jdGlvbiBvcmllbnQoKXt2YXIgdHVwbGU9dGhpcy50dXBsZTtyZXR1cm4gdGVzdChcIl1cbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIGlmKGkgPiAwKSB7XG4gICAgICBjb2RlLnB1c2goXCIsXCIpXG4gICAgfVxuICAgIGNvZGUucHVzaChcInR1cGxlW1wiLCBpLCBcIl1cIilcbiAgfVxuICBjb2RlLnB1c2goXCIpfXJldHVybiBvcmllbnRcIilcbiAgdmFyIHByb2MgPSBuZXcgRnVuY3Rpb24oXCJ0ZXN0XCIsIGNvZGUuam9pbihcIlwiKSlcbiAgdmFyIHRlc3QgPSBvcmllbnRbZCsxXVxuICBpZighdGVzdCkge1xuICAgIHRlc3QgPSBvcmllbnRcbiAgfVxuICByZXR1cm4gcHJvYyh0ZXN0KVxufVxuXG52YXIgQkFLRUQgPSBbXVxuXG5mdW5jdGlvbiBUcmlhbmd1bGF0aW9uKGRpbWVuc2lvbiwgdmVydGljZXMsIHNpbXBsaWNlcykge1xuICB0aGlzLmRpbWVuc2lvbiA9IGRpbWVuc2lvblxuICB0aGlzLnZlcnRpY2VzID0gdmVydGljZXNcbiAgdGhpcy5zaW1wbGljZXMgPSBzaW1wbGljZXNcbiAgdGhpcy5pbnRlcmlvciA9IHNpbXBsaWNlcy5maWx0ZXIoZnVuY3Rpb24oYykge1xuICAgIHJldHVybiAhYy5ib3VuZGFyeVxuICB9KVxuXG4gIHRoaXMudHVwbGUgPSBuZXcgQXJyYXkoZGltZW5zaW9uKzEpXG4gIGZvcih2YXIgaT0wOyBpPD1kaW1lbnNpb247ICsraSkge1xuICAgIHRoaXMudHVwbGVbaV0gPSB0aGlzLnZlcnRpY2VzW2ldXG4gIH1cblxuICB2YXIgbyA9IEJBS0VEW2RpbWVuc2lvbl1cbiAgaWYoIW8pIHtcbiAgICBvID0gQkFLRURbZGltZW5zaW9uXSA9IGJha2VPcmllbnQoZGltZW5zaW9uKVxuICB9XG4gIHRoaXMub3JpZW50ID0gb1xufVxuXG52YXIgcHJvdG8gPSBUcmlhbmd1bGF0aW9uLnByb3RvdHlwZVxuXG4vL0RlZ2VuZXJhdGUgc2l0dWF0aW9uIHdoZXJlIHdlIGFyZSBvbiBib3VuZGFyeSwgYnV0IGNvcGxhbmFyIHRvIGZhY2VcbnByb3RvLmhhbmRsZUJvdW5kYXJ5RGVnZW5lcmFjeSA9IGZ1bmN0aW9uKGNlbGwsIHBvaW50KSB7XG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIG4gPSB0aGlzLnZlcnRpY2VzLmxlbmd0aCAtIDFcbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG5cbiAgLy9EdW1iIHNvbHV0aW9uOiBKdXN0IGRvIGRmcyBmcm9tIGJvdW5kYXJ5IGNlbGwgdW50aWwgd2UgZmluZCBhbnkgcGVhaywgb3IgdGVybWluYXRlXG4gIHZhciB0b1Zpc2l0ID0gWyBjZWxsIF1cbiAgY2VsbC5sYXN0VmlzaXRlZCA9IC1uXG4gIHdoaWxlKHRvVmlzaXQubGVuZ3RoID4gMCkge1xuICAgIGNlbGwgPSB0b1Zpc2l0LnBvcCgpXG4gICAgdmFyIGNlbGxWZXJ0cyA9IGNlbGwudmVydGljZXNcbiAgICB2YXIgY2VsbEFkaiA9IGNlbGwuYWRqYWNlbnRcbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICB2YXIgbmVpZ2hib3IgPSBjZWxsQWRqW2ldXG4gICAgICBpZighbmVpZ2hib3IuYm91bmRhcnkgfHwgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPD0gLW4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHZhciBudiA9IG5laWdoYm9yLnZlcnRpY2VzXG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIHZhciB2diA9IG52W2pdXG4gICAgICAgIGlmKHZ2IDwgMCkge1xuICAgICAgICAgIHR1cGxlW2pdID0gcG9pbnRcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICB0dXBsZVtqXSA9IHZlcnRzW3Z2XVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICB2YXIgbyA9IHRoaXMub3JpZW50KClcbiAgICAgIGlmKG8gPiAwKSB7XG4gICAgICAgIHJldHVybiBuZWlnaGJvclxuICAgICAgfVxuICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSAtblxuICAgICAgaWYobyA9PT0gMCkge1xuICAgICAgICB0b1Zpc2l0LnB1c2gobmVpZ2hib3IpXG4gICAgICB9XG4gICAgfVxuICB9XG4gIHJldHVybiBudWxsXG59XG5cbnByb3RvLndhbGsgPSBmdW5jdGlvbihwb2ludCwgcmFuZG9tKSB7XG4gIC8vQWxpYXMgbG9jYWwgcHJvcGVydGllc1xuICB2YXIgbiA9IHRoaXMudmVydGljZXMubGVuZ3RoIC0gMVxuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuXG4gIC8vQ29tcHV0ZSBpbml0aWFsIGp1bXAgY2VsbFxuICB2YXIgaW5pdEluZGV4ID0gcmFuZG9tID8gKHRoaXMuaW50ZXJpb3IubGVuZ3RoICogTWF0aC5yYW5kb20oKSl8MCA6ICh0aGlzLmludGVyaW9yLmxlbmd0aC0xKVxuICB2YXIgY2VsbCA9IHRoaXMuaW50ZXJpb3JbIGluaXRJbmRleCBdXG5cbiAgLy9TdGFydCB3YWxraW5nXG5vdXRlckxvb3A6XG4gIHdoaWxlKCFjZWxsLmJvdW5kYXJ5KSB7XG4gICAgdmFyIGNlbGxWZXJ0cyA9IGNlbGwudmVydGljZXNcbiAgICB2YXIgY2VsbEFkaiA9IGNlbGwuYWRqYWNlbnRcblxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIHR1cGxlW2ldID0gdmVydHNbY2VsbFZlcnRzW2ldXVxuICAgIH1cbiAgICBjZWxsLmxhc3RWaXNpdGVkID0gblxuXG4gICAgLy9GaW5kIGZhcnRoZXN0IGFkamFjZW50IGNlbGxcbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICB2YXIgbmVpZ2hib3IgPSBjZWxsQWRqW2ldXG4gICAgICBpZihuZWlnaGJvci5sYXN0VmlzaXRlZCA+PSBuKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICB2YXIgcHJldiA9IHR1cGxlW2ldXG4gICAgICB0dXBsZVtpXSA9IHBvaW50XG4gICAgICB2YXIgbyA9IHRoaXMub3JpZW50KClcbiAgICAgIHR1cGxlW2ldID0gcHJldlxuICAgICAgaWYobyA8IDApIHtcbiAgICAgICAgY2VsbCA9IG5laWdoYm9yXG4gICAgICAgIGNvbnRpbnVlIG91dGVyTG9vcFxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgaWYoIW5laWdoYm9yLmJvdW5kYXJ5KSB7XG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSBuXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSAtblxuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICAgIHJldHVyblxuICB9XG5cbiAgcmV0dXJuIGNlbGxcbn1cblxucHJvdG8uYWRkUGVha3MgPSBmdW5jdGlvbihwb2ludCwgY2VsbCkge1xuICB2YXIgbiA9IHRoaXMudmVydGljZXMubGVuZ3RoIC0gMVxuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuICB2YXIgaW50ZXJpb3IgPSB0aGlzLmludGVyaW9yXG4gIHZhciBzaW1wbGljZXMgPSB0aGlzLnNpbXBsaWNlc1xuXG4gIC8vV2Fsa2luZyBmaW5pc2hlZCBhdCBib3VuZGFyeSwgdGltZSB0byBhZGQgcGVha3NcbiAgdmFyIHRvdmlzaXQgPSBbIGNlbGwgXVxuXG4gIC8vU3RyZXRjaCBpbml0aWFsIGJvdW5kYXJ5IGNlbGwgaW50byBhIHBlYWtcbiAgY2VsbC5sYXN0VmlzaXRlZCA9IG5cbiAgY2VsbC52ZXJ0aWNlc1tjZWxsLnZlcnRpY2VzLmluZGV4T2YoLTEpXSA9IG5cbiAgY2VsbC5ib3VuZGFyeSA9IGZhbHNlXG4gIGludGVyaW9yLnB1c2goY2VsbClcblxuICAvL1JlY29yZCBhIGxpc3Qgb2YgYWxsIG5ldyBib3VuZGFyaWVzIGNyZWF0ZWQgYnkgYWRkZWQgcGVha3Mgc28gd2UgY2FuIGdsdWUgdGhlbSB0b2dldGhlciB3aGVuIHdlIGFyZSBhbGwgZG9uZVxuICB2YXIgZ2x1ZUZhY2V0cyA9IFtdXG5cbiAgLy9EbyBhIHRyYXZlcnNhbCBvZiB0aGUgYm91bmRhcnkgd2Fsa2luZyBvdXR3YXJkIGZyb20gc3RhcnRpbmcgcGVha1xuICB3aGlsZSh0b3Zpc2l0Lmxlbmd0aCA+IDApIHtcbiAgICAvL1BvcCBvZmYgcGVhayBhbmQgd2FsayBvdmVyIGFkamFjZW50IGNlbGxzXG4gICAgdmFyIGNlbGwgPSB0b3Zpc2l0LnBvcCgpXG4gICAgdmFyIGNlbGxWZXJ0cyA9IGNlbGwudmVydGljZXNcbiAgICB2YXIgY2VsbEFkaiA9IGNlbGwuYWRqYWNlbnRcbiAgICB2YXIgaW5kZXhPZk4gPSBjZWxsVmVydHMuaW5kZXhPZihuKVxuICAgIGlmKGluZGV4T2ZOIDwgMCkge1xuICAgICAgY29udGludWVcbiAgICB9XG5cbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICBpZihpID09PSBpbmRleE9mTikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuXG4gICAgICAvL0ZvciBlYWNoIGJvdW5kYXJ5IG5laWdoYm9yIG9mIHRoZSBjZWxsXG4gICAgICB2YXIgbmVpZ2hib3IgPSBjZWxsQWRqW2ldXG4gICAgICBpZighbmVpZ2hib3IuYm91bmRhcnkgfHwgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPj0gbikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuXG4gICAgICB2YXIgbnYgPSBuZWlnaGJvci52ZXJ0aWNlc1xuXG4gICAgICAvL1Rlc3QgaWYgbmVpZ2hib3IgaXMgYSBwZWFrXG4gICAgICBpZihuZWlnaGJvci5sYXN0VmlzaXRlZCAhPT0gLW4pIHsgICAgICBcbiAgICAgICAgLy9Db21wdXRlIG9yaWVudGF0aW9uIG9mIHAgcmVsYXRpdmUgdG8gZWFjaCBib3VuZGFyeSBwZWFrXG4gICAgICAgIHZhciBpbmRleE9mTmVnMSA9IDBcbiAgICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICAgIGlmKG52W2pdIDwgMCkge1xuICAgICAgICAgICAgaW5kZXhPZk5lZzEgPSBqXG4gICAgICAgICAgICB0dXBsZVtqXSA9IHBvaW50XG4gICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHR1cGxlW2pdID0gdmVydHNbbnZbal1dXG4gICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHZhciBvID0gdGhpcy5vcmllbnQoKVxuXG4gICAgICAgIC8vVGVzdCBpZiBuZWlnaGJvciBjZWxsIGlzIGFsc28gYSBwZWFrXG4gICAgICAgIGlmKG8gPiAwKSB7XG4gICAgICAgICAgbnZbaW5kZXhPZk5lZzFdID0gblxuICAgICAgICAgIG5laWdoYm9yLmJvdW5kYXJ5ID0gZmFsc2VcbiAgICAgICAgICBpbnRlcmlvci5wdXNoKG5laWdoYm9yKVxuICAgICAgICAgIHRvdmlzaXQucHVzaChuZWlnaGJvcilcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IG5cbiAgICAgICAgICBjb250aW51ZVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gLW5cbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgICB2YXIgbmEgPSBuZWlnaGJvci5hZGphY2VudFxuXG4gICAgICAvL090aGVyd2lzZSwgcmVwbGFjZSBuZWlnaGJvciB3aXRoIG5ldyBmYWNlXG4gICAgICB2YXIgdnZlcnRzID0gY2VsbFZlcnRzLnNsaWNlKClcbiAgICAgIHZhciB2YWRqID0gY2VsbEFkai5zbGljZSgpXG4gICAgICB2YXIgbmNlbGwgPSBuZXcgU2ltcGxleCh2dmVydHMsIHZhZGosIHRydWUpXG4gICAgICBzaW1wbGljZXMucHVzaChuY2VsbClcblxuICAgICAgLy9Db25uZWN0IHRvIG5laWdoYm9yXG4gICAgICB2YXIgb3Bwb3NpdGUgPSBuYS5pbmRleE9mKGNlbGwpXG4gICAgICBpZihvcHBvc2l0ZSA8IDApIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIG5hW29wcG9zaXRlXSA9IG5jZWxsXG4gICAgICB2YWRqW2luZGV4T2ZOXSA9IG5laWdoYm9yXG5cbiAgICAgIC8vQ29ubmVjdCB0byBjZWxsXG4gICAgICB2dmVydHNbaV0gPSAtMVxuICAgICAgdmFkaltpXSA9IGNlbGxcbiAgICAgIGNlbGxBZGpbaV0gPSBuY2VsbFxuXG4gICAgICAvL0ZsaXAgZmFjZXRcbiAgICAgIG5jZWxsLmZsaXAoKVxuXG4gICAgICAvL0FkZCB0byBnbHVlIGxpc3RcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgdmFyIHV1ID0gdnZlcnRzW2pdXG4gICAgICAgIGlmKHV1IDwgMCB8fCB1dSA9PT0gbikge1xuICAgICAgICAgIGNvbnRpbnVlXG4gICAgICAgIH1cbiAgICAgICAgdmFyIG5mYWNlID0gbmV3IEFycmF5KGQtMSlcbiAgICAgICAgdmFyIG5wdHIgPSAwXG4gICAgICAgIGZvcih2YXIgaz0wOyBrPD1kOyArK2spIHtcbiAgICAgICAgICB2YXIgdnYgPSB2dmVydHNba11cbiAgICAgICAgICBpZih2diA8IDAgfHwgayA9PT0gaikge1xuICAgICAgICAgICAgY29udGludWVcbiAgICAgICAgICB9XG4gICAgICAgICAgbmZhY2VbbnB0cisrXSA9IHZ2XG4gICAgICAgIH1cbiAgICAgICAgZ2x1ZUZhY2V0cy5wdXNoKG5ldyBHbHVlRmFjZXQobmZhY2UsIG5jZWxsLCBqKSlcbiAgICAgIH1cbiAgICB9XG4gIH1cblxuICAvL0dsdWUgYm91bmRhcnkgZmFjZXRzIHRvZ2V0aGVyXG4gIGdsdWVGYWNldHMuc29ydChjb21wYXJlR2x1ZSlcblxuICBmb3IodmFyIGk9MDsgaSsxPGdsdWVGYWNldHMubGVuZ3RoOyBpKz0yKSB7XG4gICAgdmFyIGEgPSBnbHVlRmFjZXRzW2ldXG4gICAgdmFyIGIgPSBnbHVlRmFjZXRzW2krMV1cbiAgICB2YXIgYWkgPSBhLmluZGV4XG4gICAgdmFyIGJpID0gYi5pbmRleFxuICAgIGlmKGFpIDwgMCB8fCBiaSA8IDApIHtcbiAgICAgIGNvbnRpbnVlXG4gICAgfVxuICAgIGEuY2VsbC5hZGphY2VudFthLmluZGV4XSA9IGIuY2VsbFxuICAgIGIuY2VsbC5hZGphY2VudFtiLmluZGV4XSA9IGEuY2VsbFxuICB9XG59XG5cbnByb3RvLmluc2VydCA9IGZ1bmN0aW9uKHBvaW50LCByYW5kb20pIHtcbiAgLy9BZGQgcG9pbnRcbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuICB2ZXJ0cy5wdXNoKHBvaW50KVxuXG4gIHZhciBjZWxsID0gdGhpcy53YWxrKHBvaW50LCByYW5kb20pXG4gIGlmKCFjZWxsKSB7XG4gICAgcmV0dXJuXG4gIH1cblxuICAvL0FsaWFzIGxvY2FsIHByb3BlcnRpZXNcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG5cbiAgLy9EZWdlbmVyYXRlIGNhc2U6IElmIHBvaW50IGlzIGNvcGxhbmFyIHRvIGNlbGwsIHRoZW4gd2FsayB1bnRpbCB3ZSBmaW5kIGEgbm9uLWRlZ2VuZXJhdGUgYm91bmRhcnlcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB2diA9IGNlbGwudmVydGljZXNbaV1cbiAgICBpZih2diA8IDApIHtcbiAgICAgIHR1cGxlW2ldID0gcG9pbnRcbiAgICB9IGVsc2Uge1xuICAgICAgdHVwbGVbaV0gPSB2ZXJ0c1t2dl1cbiAgICB9XG4gIH1cbiAgdmFyIG8gPSB0aGlzLm9yaWVudCh0dXBsZSlcbiAgaWYobyA8IDApIHtcbiAgICByZXR1cm5cbiAgfSBlbHNlIGlmKG8gPT09IDApIHtcbiAgICBjZWxsID0gdGhpcy5oYW5kbGVCb3VuZGFyeURlZ2VuZXJhY3koY2VsbCwgcG9pbnQpXG4gICAgaWYoIWNlbGwpIHtcbiAgICAgIHJldHVyblxuICAgIH1cbiAgfVxuXG4gIC8vQWRkIHBlYWtzXG4gIHRoaXMuYWRkUGVha3MocG9pbnQsIGNlbGwpXG59XG5cbi8vRXh0cmFjdCBhbGwgYm91bmRhcnkgY2VsbHNcbnByb3RvLmJvdW5kYXJ5ID0gZnVuY3Rpb24oKSB7XG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIGJvdW5kYXJ5ID0gW11cbiAgdmFyIGNlbGxzID0gdGhpcy5zaW1wbGljZXNcbiAgdmFyIG5jID0gY2VsbHMubGVuZ3RoXG4gIGZvcih2YXIgaT0wOyBpPG5jOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgaWYoYy5ib3VuZGFyeSkge1xuICAgICAgdmFyIGJjZWxsID0gbmV3IEFycmF5KGQpXG4gICAgICB2YXIgY3YgPSBjLnZlcnRpY2VzXG4gICAgICB2YXIgcHRyID0gMFxuICAgICAgdmFyIHBhcml0eSA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgaWYoY3Zbal0gPj0gMCkge1xuICAgICAgICAgIGJjZWxsW3B0cisrXSA9IGN2W2pdXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgcGFyaXR5ID0gaiYxXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIGlmKHBhcml0eSA9PT0gKGQmMSkpIHtcbiAgICAgICAgdmFyIHQgPSBiY2VsbFswXVxuICAgICAgICBiY2VsbFswXSA9IGJjZWxsWzFdXG4gICAgICAgIGJjZWxsWzFdID0gdFxuICAgICAgfVxuICAgICAgYm91bmRhcnkucHVzaChiY2VsbClcbiAgICB9XG4gIH1cbiAgcmV0dXJuIGJvdW5kYXJ5XG59XG5cbmZ1bmN0aW9uIGluY3JlbWVudGFsQ29udmV4SHVsbChwb2ludHMsIHJhbmRvbVNlYXJjaCkge1xuICB2YXIgbiA9IHBvaW50cy5sZW5ndGhcbiAgaWYobiA9PT0gMCkge1xuICAgIHRocm93IG5ldyBFcnJvcihcIk11c3QgaGF2ZSBhdCBsZWFzdCBkKzEgcG9pbnRzXCIpXG4gIH1cbiAgdmFyIGQgPSBwb2ludHNbMF0ubGVuZ3RoXG4gIGlmKG4gPD0gZCkge1xuICAgIHRocm93IG5ldyBFcnJvcihcIk11c3QgaW5wdXQgYXQgbGVhc3QgZCsxIHBvaW50c1wiKVxuICB9XG5cbiAgLy9GSVhNRTogVGhpcyBjb3VsZCBiZSBkZWdlbmVyYXRlLCBidXQgbmVlZCB0byBzZWxlY3QgZCsxIG5vbi1jb3BsYW5hciBwb2ludHMgdG8gYm9vdHN0cmFwIHByb2Nlc3NcbiAgdmFyIGluaXRpYWxTaW1wbGV4ID0gcG9pbnRzLnNsaWNlKDAsIGQrMSlcblxuICAvL01ha2Ugc3VyZSBpbml0aWFsIHNpbXBsZXggaXMgcG9zaXRpdmVseSBvcmllbnRlZFxuICB2YXIgbyA9IG9yaWVudC5hcHBseSh2b2lkIDAsIGluaXRpYWxTaW1wbGV4KVxuICBpZihvID09PSAwKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKFwiSW5wdXQgbm90IGluIGdlbmVyYWwgcG9zaXRpb25cIilcbiAgfVxuICB2YXIgaW5pdGlhbENvb3JkcyA9IG5ldyBBcnJheShkKzEpXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICBpbml0aWFsQ29vcmRzW2ldID0gaVxuICB9XG4gIGlmKG8gPCAwKSB7XG4gICAgaW5pdGlhbENvb3Jkc1swXSA9IDFcbiAgICBpbml0aWFsQ29vcmRzWzFdID0gMFxuICB9XG5cbiAgLy9DcmVhdGUgaW5pdGlhbCB0b3BvbG9naWNhbCBpbmRleCwgZ2x1ZSBwb2ludGVycyB0b2dldGhlciAoa2luZCBvZiBtZXNzeSlcbiAgdmFyIGluaXRpYWxDZWxsID0gbmV3IFNpbXBsZXgoaW5pdGlhbENvb3JkcywgbmV3IEFycmF5KGQrMSksIGZhbHNlKVxuICB2YXIgYm91bmRhcnkgPSBpbml0aWFsQ2VsbC5hZGphY2VudFxuICB2YXIgbGlzdCA9IG5ldyBBcnJheShkKzIpXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgdmVydHMgPSBpbml0aWFsQ29vcmRzLnNsaWNlKClcbiAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICBpZihqID09PSBpKSB7XG4gICAgICAgIHZlcnRzW2pdID0gLTFcbiAgICAgIH1cbiAgICB9XG4gICAgdmFyIHQgPSB2ZXJ0c1swXVxuICAgIHZlcnRzWzBdID0gdmVydHNbMV1cbiAgICB2ZXJ0c1sxXSA9IHRcbiAgICB2YXIgY2VsbCA9IG5ldyBTaW1wbGV4KHZlcnRzLCBuZXcgQXJyYXkoZCsxKSwgdHJ1ZSlcbiAgICBib3VuZGFyeVtpXSA9IGNlbGxcbiAgICBsaXN0W2ldID0gY2VsbFxuICB9XG4gIGxpc3RbZCsxXSA9IGluaXRpYWxDZWxsXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgdmVydHMgPSBib3VuZGFyeVtpXS52ZXJ0aWNlc1xuICAgIHZhciBhZGogPSBib3VuZGFyeVtpXS5hZGphY2VudFxuICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgIHZhciB2ID0gdmVydHNbal1cbiAgICAgIGlmKHYgPCAwKSB7XG4gICAgICAgIGFkaltqXSA9IGluaXRpYWxDZWxsXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBmb3IodmFyIGs9MDsgazw9ZDsgKytrKSB7XG4gICAgICAgIGlmKGJvdW5kYXJ5W2tdLnZlcnRpY2VzLmluZGV4T2YodikgPCAwKSB7XG4gICAgICAgICAgYWRqW2pdID0gYm91bmRhcnlba11cbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgfVxuXG4gIC8vSW5pdGlhbGl6ZSB0cmlhbmdsZXNcbiAgdmFyIHRyaWFuZ2xlcyA9IG5ldyBUcmlhbmd1bGF0aW9uKGQsIGluaXRpYWxTaW1wbGV4LCBsaXN0KVxuXG4gIC8vSW5zZXJ0IHJlbWFpbmluZyBwb2ludHNcbiAgdmFyIHVzZVJhbmRvbSA9ICEhcmFuZG9tU2VhcmNoXG4gIGZvcih2YXIgaT1kKzE7IGk8bjsgKytpKSB7XG4gICAgdHJpYW5nbGVzLmluc2VydChwb2ludHNbaV0sIHVzZVJhbmRvbSlcbiAgfVxuICBcbiAgLy9FeHRyYWN0IGJvdW5kYXJ5IGNlbGxzXG4gIHJldHVybiB0cmlhbmdsZXMuYm91bmRhcnkoKVxufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gZmFzdFR3b1N1bVxuXG5mdW5jdGlvbiBmYXN0VHdvU3VtKGEsIGIsIHJlc3VsdCkge1xuXHR2YXIgeCA9IGEgKyBiXG5cdHZhciBidiA9IHggLSBhXG5cdHZhciBhdiA9IHggLSBidlxuXHR2YXIgYnIgPSBiIC0gYnZcblx0dmFyIGFyID0gYSAtIGF2XG5cdGlmKHJlc3VsdCkge1xuXHRcdHJlc3VsdFswXSA9IGFyICsgYnJcblx0XHRyZXN1bHRbMV0gPSB4XG5cdFx0cmV0dXJuIHJlc3VsdFxuXHR9XG5cdHJldHVybiBbYXIrYnIsIHhdXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxudmFyIHR3b1Byb2R1Y3QgPSByZXF1aXJlKFwidHdvLXByb2R1Y3RcIilcbnZhciB0d29TdW0gPSByZXF1aXJlKFwidHdvLXN1bVwiKVxuXG5tb2R1bGUuZXhwb3J0cyA9IHNjYWxlTGluZWFyRXhwYW5zaW9uXG5cbmZ1bmN0aW9uIHNjYWxlTGluZWFyRXhwYW5zaW9uKGUsIHNjYWxlKSB7XG4gIHZhciBuID0gZS5sZW5ndGhcbiAgaWYobiA9PT0gMSkge1xuICAgIHZhciB0cyA9IHR3b1Byb2R1Y3QoZVswXSwgc2NhbGUpXG4gICAgaWYodHNbMF0pIHtcbiAgICAgIHJldHVybiB0c1xuICAgIH1cbiAgICByZXR1cm4gWyB0c1sxXSBdXG4gIH1cbiAgdmFyIGcgPSBuZXcgQXJyYXkoMiAqIG4pXG4gIHZhciBxID0gWzAuMSwgMC4xXVxuICB2YXIgdCA9IFswLjEsIDAuMV1cbiAgdmFyIGNvdW50ID0gMFxuICB0d29Qcm9kdWN0KGVbMF0sIHNjYWxlLCBxKVxuICBpZihxWzBdKSB7XG4gICAgZ1tjb3VudCsrXSA9IHFbMF1cbiAgfVxuICBmb3IodmFyIGk9MTsgaTxuOyArK2kpIHtcbiAgICB0d29Qcm9kdWN0KGVbaV0sIHNjYWxlLCB0KVxuICAgIHZhciBwcSA9IHFbMV1cbiAgICB0d29TdW0ocHEsIHRbMF0sIHEpXG4gICAgaWYocVswXSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHFbMF1cbiAgICB9XG4gICAgdmFyIGEgPSB0WzFdXG4gICAgdmFyIGIgPSBxWzFdXG4gICAgdmFyIHggPSBhICsgYlxuICAgIHZhciBidiA9IHggLSBhXG4gICAgdmFyIHkgPSBiIC0gYnZcbiAgICBxWzFdID0geFxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICB9XG4gIGlmKHFbMV0pIHtcbiAgICBnW2NvdW50KytdID0gcVsxXVxuICB9XG4gIGlmKGNvdW50ID09PSAwKSB7XG4gICAgZ1tjb3VudCsrXSA9IDAuMFxuICB9XG4gIGcubGVuZ3RoID0gY291bnRcbiAgcmV0dXJuIGdcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IHJvYnVzdFN1YnRyYWN0XG5cbi8vRWFzeSBjYXNlOiBBZGQgdHdvIHNjYWxhcnNcbmZ1bmN0aW9uIHNjYWxhclNjYWxhcihhLCBiKSB7XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIGF2ID0geCAtIGJ2XG4gIHZhciBiciA9IGIgLSBidlxuICB2YXIgYXIgPSBhIC0gYXZcbiAgdmFyIHkgPSBhciArIGJyXG4gIGlmKHkpIHtcbiAgICByZXR1cm4gW3ksIHhdXG4gIH1cbiAgcmV0dXJuIFt4XVxufVxuXG5mdW5jdGlvbiByb2J1c3RTdWJ0cmFjdChlLCBmKSB7XG4gIHZhciBuZSA9IGUubGVuZ3RofDBcbiAgdmFyIG5mID0gZi5sZW5ndGh8MFxuICBpZihuZSA9PT0gMSAmJiBuZiA9PT0gMSkge1xuICAgIHJldHVybiBzY2FsYXJTY2FsYXIoZVswXSwgLWZbMF0pXG4gIH1cbiAgdmFyIG4gPSBuZSArIG5mXG4gIHZhciBnID0gbmV3IEFycmF5KG4pXG4gIHZhciBjb3VudCA9IDBcbiAgdmFyIGVwdHIgPSAwXG4gIHZhciBmcHRyID0gMFxuICB2YXIgYWJzID0gTWF0aC5hYnNcbiAgdmFyIGVpID0gZVtlcHRyXVxuICB2YXIgZWEgPSBhYnMoZWkpXG4gIHZhciBmaSA9IC1mW2ZwdHJdXG4gIHZhciBmYSA9IGFicyhmaSlcbiAgdmFyIGEsIGJcbiAgaWYoZWEgPCBmYSkge1xuICAgIGIgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYiA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgaWYoKGVwdHIgPCBuZSAmJiBlYSA8IGZhKSB8fCAoZnB0ciA+PSBuZikpIHtcbiAgICBhID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGEgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSAtZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIHkgPSBiIC0gYnZcbiAgdmFyIHEwID0geVxuICB2YXIgcTEgPSB4XG4gIHZhciBfeCwgX2J2LCBfYXYsIF9iciwgX2FyXG4gIHdoaWxlKGVwdHIgPCBuZSAmJiBmcHRyIDwgbmYpIHtcbiAgICBpZihlYSA8IGZhKSB7XG4gICAgICBhID0gZWlcbiAgICAgIGVwdHIgKz0gMVxuICAgICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgICBlYSA9IGFicyhlaSlcbiAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgYSA9IGZpXG4gICAgICBmcHRyICs9IDFcbiAgICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgICAgIGZhID0gYWJzKGZpKVxuICAgICAgfVxuICAgIH1cbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gIH1cbiAgd2hpbGUoZXB0ciA8IG5lKSB7XG4gICAgYSA9IGVpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgfVxuICB9XG4gIHdoaWxlKGZwdHIgPCBuZikge1xuICAgIGEgPSBmaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9IFxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgfVxuICB9XG4gIGlmKHEwKSB7XG4gICAgZ1tjb3VudCsrXSA9IHEwXG4gIH1cbiAgaWYocTEpIHtcbiAgICBnW2NvdW50KytdID0gcTFcbiAgfVxuICBpZighY291bnQpIHtcbiAgICBnW2NvdW50KytdID0gMC4wICBcbiAgfVxuICBnLmxlbmd0aCA9IGNvdW50XG4gIHJldHVybiBnXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSBsaW5lYXJFeHBhbnNpb25TdW1cblxuLy9FYXN5IGNhc2U6IEFkZCB0d28gc2NhbGFyc1xuZnVuY3Rpb24gc2NhbGFyU2NhbGFyKGEsIGIpIHtcbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgYXYgPSB4IC0gYnZcbiAgdmFyIGJyID0gYiAtIGJ2XG4gIHZhciBhciA9IGEgLSBhdlxuICB2YXIgeSA9IGFyICsgYnJcbiAgaWYoeSkge1xuICAgIHJldHVybiBbeSwgeF1cbiAgfVxuICByZXR1cm4gW3hdXG59XG5cbmZ1bmN0aW9uIGxpbmVhckV4cGFuc2lvblN1bShlLCBmKSB7XG4gIHZhciBuZSA9IGUubGVuZ3RofDBcbiAgdmFyIG5mID0gZi5sZW5ndGh8MFxuICBpZihuZSA9PT0gMSAmJiBuZiA9PT0gMSkge1xuICAgIHJldHVybiBzY2FsYXJTY2FsYXIoZVswXSwgZlswXSlcbiAgfVxuICB2YXIgbiA9IG5lICsgbmZcbiAgdmFyIGcgPSBuZXcgQXJyYXkobilcbiAgdmFyIGNvdW50ID0gMFxuICB2YXIgZXB0ciA9IDBcbiAgdmFyIGZwdHIgPSAwXG4gIHZhciBhYnMgPSBNYXRoLmFic1xuICB2YXIgZWkgPSBlW2VwdHJdXG4gIHZhciBlYSA9IGFicyhlaSlcbiAgdmFyIGZpID0gZltmcHRyXVxuICB2YXIgZmEgPSBhYnMoZmkpXG4gIHZhciBhLCBiXG4gIGlmKGVhIDwgZmEpIHtcbiAgICBiID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGIgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgaWYoKGVwdHIgPCBuZSAmJiBlYSA8IGZhKSB8fCAoZnB0ciA+PSBuZikpIHtcbiAgICBhID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGEgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgeSA9IGIgLSBidlxuICB2YXIgcTAgPSB5XG4gIHZhciBxMSA9IHhcbiAgdmFyIF94LCBfYnYsIF9hdiwgX2JyLCBfYXJcbiAgd2hpbGUoZXB0ciA8IG5lICYmIGZwdHIgPCBuZikge1xuICAgIGlmKGVhIDwgZmEpIHtcbiAgICAgIGEgPSBlaVxuICAgICAgZXB0ciArPSAxXG4gICAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICAgIGVhID0gYWJzKGVpKVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICBhID0gZmlcbiAgICAgIGZwdHIgKz0gMVxuICAgICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICAgIGZpID0gZltmcHRyXVxuICAgICAgICBmYSA9IGFicyhmaSlcbiAgICAgIH1cbiAgICB9XG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICB9XG4gIHdoaWxlKGVwdHIgPCBuZSkge1xuICAgIGEgPSBlaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgIH1cbiAgfVxuICB3aGlsZShmcHRyIDwgbmYpIHtcbiAgICBhID0gZmlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfSBcbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgfVxuICB9XG4gIGlmKHEwKSB7XG4gICAgZ1tjb3VudCsrXSA9IHEwXG4gIH1cbiAgaWYocTEpIHtcbiAgICBnW2NvdW50KytdID0gcTFcbiAgfVxuICBpZighY291bnQpIHtcbiAgICBnW2NvdW50KytdID0gMC4wICBcbiAgfVxuICBnLmxlbmd0aCA9IGNvdW50XG4gIHJldHVybiBnXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSB0d29Qcm9kdWN0XG5cbnZhciBTUExJVFRFUiA9ICsoTWF0aC5wb3coMiwgMjcpICsgMS4wKVxuXG5mdW5jdGlvbiB0d29Qcm9kdWN0KGEsIGIsIHJlc3VsdCkge1xuICB2YXIgeCA9IGEgKiBiXG5cbiAgdmFyIGMgPSBTUExJVFRFUiAqIGFcbiAgdmFyIGFiaWcgPSBjIC0gYVxuICB2YXIgYWhpID0gYyAtIGFiaWdcbiAgdmFyIGFsbyA9IGEgLSBhaGlcblxuICB2YXIgZCA9IFNQTElUVEVSICogYlxuICB2YXIgYmJpZyA9IGQgLSBiXG4gIHZhciBiaGkgPSBkIC0gYmJpZ1xuICB2YXIgYmxvID0gYiAtIGJoaVxuXG4gIHZhciBlcnIxID0geCAtIChhaGkgKiBiaGkpXG4gIHZhciBlcnIyID0gZXJyMSAtIChhbG8gKiBiaGkpXG4gIHZhciBlcnIzID0gZXJyMiAtIChhaGkgKiBibG8pXG5cbiAgdmFyIHkgPSBhbG8gKiBibG8gLSBlcnIzXG5cbiAgaWYocmVzdWx0KSB7XG4gICAgcmVzdWx0WzBdID0geVxuICAgIHJlc3VsdFsxXSA9IHhcbiAgICByZXR1cm4gcmVzdWx0XG4gIH1cblxuICByZXR1cm4gWyB5LCB4IF1cbn0iLCJcInVzZSBzdHJpY3RcIlxuXG52YXIgdHdvUHJvZHVjdCA9IHJlcXVpcmUoXCJ0d28tcHJvZHVjdFwiKVxudmFyIHJvYnVzdFN1bSA9IHJlcXVpcmUoXCJyb2J1c3Qtc3VtXCIpXG52YXIgcm9idXN0U2NhbGUgPSByZXF1aXJlKFwicm9idXN0LXNjYWxlXCIpXG52YXIgcm9idXN0U3VidHJhY3QgPSByZXF1aXJlKFwicm9idXN0LXN1YnRyYWN0XCIpXG5cbnZhciBOVU1fRVhQQU5EID0gNVxuXG52YXIgRVBTSUxPTiAgICAgPSAxLjExMDIyMzAyNDYyNTE1NjVlLTE2XG52YXIgRVJSQk9VTkQzICAgPSAoMy4wICsgMTYuMCAqIEVQU0lMT04pICogRVBTSUxPTlxudmFyIEVSUkJPVU5ENCAgID0gKDcuMCArIDU2LjAgKiBFUFNJTE9OKSAqIEVQU0lMT05cblxuZnVuY3Rpb24gY29mYWN0b3IobSwgYykge1xuICB2YXIgcmVzdWx0ID0gbmV3IEFycmF5KG0ubGVuZ3RoLTEpXG4gIGZvcih2YXIgaT0xOyBpPG0ubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgciA9IHJlc3VsdFtpLTFdID0gbmV3IEFycmF5KG0ubGVuZ3RoLTEpXG4gICAgZm9yKHZhciBqPTAsaz0wOyBqPG0ubGVuZ3RoOyArK2opIHtcbiAgICAgIGlmKGogPT09IGMpIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHJbaysrXSA9IG1baV1bal1cbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJlc3VsdFxufVxuXG5mdW5jdGlvbiBtYXRyaXgobikge1xuICB2YXIgcmVzdWx0ID0gbmV3IEFycmF5KG4pXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIHJlc3VsdFtpXSA9IG5ldyBBcnJheShuKVxuICAgIGZvcih2YXIgaj0wOyBqPG47ICsraikge1xuICAgICAgcmVzdWx0W2ldW2pdID0gW1wibVwiLCBqLCBcIltcIiwgKG4taS0xKSwgXCJdXCJdLmpvaW4oXCJcIilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJlc3VsdFxufVxuXG5mdW5jdGlvbiBzaWduKG4pIHtcbiAgaWYobiAmIDEpIHtcbiAgICByZXR1cm4gXCItXCJcbiAgfVxuICByZXR1cm4gXCJcIlxufVxuXG5mdW5jdGlvbiBnZW5lcmF0ZVN1bShleHByKSB7XG4gIGlmKGV4cHIubGVuZ3RoID09PSAxKSB7XG4gICAgcmV0dXJuIGV4cHJbMF1cbiAgfSBlbHNlIGlmKGV4cHIubGVuZ3RoID09PSAyKSB7XG4gICAgcmV0dXJuIFtcInN1bShcIiwgZXhwclswXSwgXCIsXCIsIGV4cHJbMV0sIFwiKVwiXS5qb2luKFwiXCIpXG4gIH0gZWxzZSB7XG4gICAgdmFyIG0gPSBleHByLmxlbmd0aD4+MVxuICAgIHJldHVybiBbXCJzdW0oXCIsIGdlbmVyYXRlU3VtKGV4cHIuc2xpY2UoMCwgbSkpLCBcIixcIiwgZ2VuZXJhdGVTdW0oZXhwci5zbGljZShtKSksIFwiKVwiXS5qb2luKFwiXCIpXG4gIH1cbn1cblxuZnVuY3Rpb24gZGV0ZXJtaW5hbnQobSkge1xuICBpZihtLmxlbmd0aCA9PT0gMikge1xuICAgIHJldHVybiBbW1wic3VtKHByb2QoXCIsIG1bMF1bMF0sIFwiLFwiLCBtWzFdWzFdLCBcIikscHJvZCgtXCIsIG1bMF1bMV0sIFwiLFwiLCBtWzFdWzBdLCBcIikpXCJdLmpvaW4oXCJcIildXG4gIH0gZWxzZSB7XG4gICAgdmFyIGV4cHIgPSBbXVxuICAgIGZvcih2YXIgaT0wOyBpPG0ubGVuZ3RoOyArK2kpIHtcbiAgICAgIGV4cHIucHVzaChbXCJzY2FsZShcIiwgZ2VuZXJhdGVTdW0oZGV0ZXJtaW5hbnQoY29mYWN0b3IobSwgaSkpKSwgXCIsXCIsIHNpZ24oaSksIG1bMF1baV0sIFwiKVwiXS5qb2luKFwiXCIpKVxuICAgIH1cbiAgICByZXR1cm4gZXhwclxuICB9XG59XG5cbmZ1bmN0aW9uIG9yaWVudGF0aW9uKG4pIHtcbiAgdmFyIHBvcyA9IFtdXG4gIHZhciBuZWcgPSBbXVxuICB2YXIgbSA9IG1hdHJpeChuKVxuICB2YXIgYXJncyA9IFtdXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIGlmKChpJjEpPT09MCkge1xuICAgICAgcG9zLnB1c2guYXBwbHkocG9zLCBkZXRlcm1pbmFudChjb2ZhY3RvcihtLCBpKSkpXG4gICAgfSBlbHNlIHtcbiAgICAgIG5lZy5wdXNoLmFwcGx5KG5lZywgZGV0ZXJtaW5hbnQoY29mYWN0b3IobSwgaSkpKVxuICAgIH1cbiAgICBhcmdzLnB1c2goXCJtXCIgKyBpKVxuICB9XG4gIHZhciBwb3NFeHByID0gZ2VuZXJhdGVTdW0ocG9zKVxuICB2YXIgbmVnRXhwciA9IGdlbmVyYXRlU3VtKG5lZylcbiAgdmFyIGZ1bmNOYW1lID0gXCJvcmllbnRhdGlvblwiICsgbiArIFwiRXhhY3RcIlxuICB2YXIgY29kZSA9IFtcImZ1bmN0aW9uIFwiLCBmdW5jTmFtZSwgXCIoXCIsIGFyZ3Muam9pbigpLCBcIil7dmFyIHA9XCIsIHBvc0V4cHIsIFwiLG49XCIsIG5lZ0V4cHIsIFwiLGQ9c3ViKHAsbik7XFxcbnJldHVybiBkW2QubGVuZ3RoLTFdO307cmV0dXJuIFwiLCBmdW5jTmFtZV0uam9pbihcIlwiKVxuICB2YXIgcHJvYyA9IG5ldyBGdW5jdGlvbihcInN1bVwiLCBcInByb2RcIiwgXCJzY2FsZVwiLCBcInN1YlwiLCBjb2RlKVxuICByZXR1cm4gcHJvYyhyb2J1c3RTdW0sIHR3b1Byb2R1Y3QsIHJvYnVzdFNjYWxlLCByb2J1c3RTdWJ0cmFjdClcbn1cblxudmFyIG9yaWVudGF0aW9uM0V4YWN0ID0gb3JpZW50YXRpb24oMylcbnZhciBvcmllbnRhdGlvbjRFeGFjdCA9IG9yaWVudGF0aW9uKDQpXG5cbnZhciBDQUNIRUQgPSBbXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMCgpIHsgcmV0dXJuIDAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24xKCkgeyByZXR1cm4gMCB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjIoYSwgYikgeyBcbiAgICByZXR1cm4gYlswXSAtIGFbMF1cbiAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24zKGEsIGIsIGMpIHtcbiAgICB2YXIgbCA9IChhWzFdIC0gY1sxXSkgKiAoYlswXSAtIGNbMF0pXG4gICAgdmFyIHIgPSAoYVswXSAtIGNbMF0pICogKGJbMV0gLSBjWzFdKVxuICAgIHZhciBkZXQgPSBsIC0gclxuICAgIHZhciBzXG4gICAgaWYobCA+IDApIHtcbiAgICAgIGlmKHIgPD0gMCkge1xuICAgICAgICByZXR1cm4gZGV0XG4gICAgICB9IGVsc2Uge1xuICAgICAgICBzID0gbCArIHJcbiAgICAgIH1cbiAgICB9IGVsc2UgaWYobCA8IDApIHtcbiAgICAgIGlmKHIgPj0gMCkge1xuICAgICAgICByZXR1cm4gZGV0XG4gICAgICB9IGVsc2Uge1xuICAgICAgICBzID0gLShsICsgcilcbiAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgcmV0dXJuIGRldFxuICAgIH1cbiAgICB2YXIgdG9sID0gRVJSQk9VTkQzICogc1xuICAgIGlmKGRldCA+PSB0b2wgfHwgZGV0IDw9IC10b2wpIHtcbiAgICAgIHJldHVybiBkZXRcbiAgICB9XG4gICAgcmV0dXJuIG9yaWVudGF0aW9uM0V4YWN0KGEsIGIsIGMpXG4gIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uNChhLGIsYyxkKSB7XG4gICAgdmFyIGFkeCA9IGFbMF0gLSBkWzBdXG4gICAgdmFyIGJkeCA9IGJbMF0gLSBkWzBdXG4gICAgdmFyIGNkeCA9IGNbMF0gLSBkWzBdXG4gICAgdmFyIGFkeSA9IGFbMV0gLSBkWzFdXG4gICAgdmFyIGJkeSA9IGJbMV0gLSBkWzFdXG4gICAgdmFyIGNkeSA9IGNbMV0gLSBkWzFdXG4gICAgdmFyIGFkeiA9IGFbMl0gLSBkWzJdXG4gICAgdmFyIGJkeiA9IGJbMl0gLSBkWzJdXG4gICAgdmFyIGNkeiA9IGNbMl0gLSBkWzJdXG4gICAgdmFyIGJkeGNkeSA9IGJkeCAqIGNkeVxuICAgIHZhciBjZHhiZHkgPSBjZHggKiBiZHlcbiAgICB2YXIgY2R4YWR5ID0gY2R4ICogYWR5XG4gICAgdmFyIGFkeGNkeSA9IGFkeCAqIGNkeVxuICAgIHZhciBhZHhiZHkgPSBhZHggKiBiZHlcbiAgICB2YXIgYmR4YWR5ID0gYmR4ICogYWR5XG4gICAgdmFyIGRldCA9IGFkeiAqIChiZHhjZHkgLSBjZHhiZHkpIFxuICAgICAgICAgICAgKyBiZHogKiAoY2R4YWR5IC0gYWR4Y2R5KVxuICAgICAgICAgICAgKyBjZHogKiAoYWR4YmR5IC0gYmR4YWR5KVxuICAgIHZhciBwZXJtYW5lbnQgPSAoTWF0aC5hYnMoYmR4Y2R5KSArIE1hdGguYWJzKGNkeGJkeSkpICogTWF0aC5hYnMoYWR6KVxuICAgICAgICAgICAgICAgICAgKyAoTWF0aC5hYnMoY2R4YWR5KSArIE1hdGguYWJzKGFkeGNkeSkpICogTWF0aC5hYnMoYmR6KVxuICAgICAgICAgICAgICAgICAgKyAoTWF0aC5hYnMoYWR4YmR5KSArIE1hdGguYWJzKGJkeGFkeSkpICogTWF0aC5hYnMoY2R6KVxuICAgIHZhciB0b2wgPSBFUlJCT1VORDQgKiBwZXJtYW5lbnRcbiAgICBpZiAoKGRldCA+IHRvbCkgfHwgKC1kZXQgPiB0b2wpKSB7XG4gICAgICByZXR1cm4gZGV0XG4gICAgfVxuICAgIHJldHVybiBvcmllbnRhdGlvbjRFeGFjdChhLGIsYyxkKVxuICB9XG5dXG5cbmZ1bmN0aW9uIHNsb3dPcmllbnQoYXJncykge1xuICB2YXIgcHJvYyA9IENBQ0hFRFthcmdzLmxlbmd0aF1cbiAgaWYoIXByb2MpIHtcbiAgICBwcm9jID0gQ0FDSEVEW2FyZ3MubGVuZ3RoXSA9IG9yaWVudGF0aW9uKGFyZ3MubGVuZ3RoKVxuICB9XG4gIHJldHVybiBwcm9jLmFwcGx5KHVuZGVmaW5lZCwgYXJncylcbn1cblxuZnVuY3Rpb24gZ2VuZXJhdGVPcmllbnRhdGlvblByb2MoKSB7XG4gIHdoaWxlKENBQ0hFRC5sZW5ndGggPD0gTlVNX0VYUEFORCkge1xuICAgIENBQ0hFRC5wdXNoKG9yaWVudGF0aW9uKENBQ0hFRC5sZW5ndGgpKVxuICB9XG4gIHZhciBhcmdzID0gW11cbiAgdmFyIHByb2NBcmdzID0gW1wic2xvd1wiXVxuICBmb3IodmFyIGk9MDsgaTw9TlVNX0VYUEFORDsgKytpKSB7XG4gICAgYXJncy5wdXNoKFwiYVwiICsgaSlcbiAgICBwcm9jQXJncy5wdXNoKFwib1wiICsgaSlcbiAgfVxuICB2YXIgY29kZSA9IFtcbiAgICBcImZ1bmN0aW9uIGdldE9yaWVudGF0aW9uKFwiLCBhcmdzLmpvaW4oKSwgXCIpe3N3aXRjaChhcmd1bWVudHMubGVuZ3RoKXtjYXNlIDA6Y2FzZSAxOnJldHVybiAwO1wiXG4gIF1cbiAgZm9yKHZhciBpPTI7IGk8PU5VTV9FWFBBTkQ7ICsraSkge1xuICAgIGNvZGUucHVzaChcImNhc2UgXCIsIGksIFwiOnJldHVybiBvXCIsIGksIFwiKFwiLCBhcmdzLnNsaWNlKDAsIGkpLmpvaW4oKSwgXCIpO1wiKVxuICB9XG4gIGNvZGUucHVzaChcIn12YXIgcz1uZXcgQXJyYXkoYXJndW1lbnRzLmxlbmd0aCk7Zm9yKHZhciBpPTA7aTxhcmd1bWVudHMubGVuZ3RoOysraSl7c1tpXT1hcmd1bWVudHNbaV19O3JldHVybiBzbG93KHMpO31yZXR1cm4gZ2V0T3JpZW50YXRpb25cIilcbiAgcHJvY0FyZ3MucHVzaChjb2RlLmpvaW4oXCJcIikpXG5cbiAgdmFyIHByb2MgPSBGdW5jdGlvbi5hcHBseSh1bmRlZmluZWQsIHByb2NBcmdzKVxuICBtb2R1bGUuZXhwb3J0cyA9IHByb2MuYXBwbHkodW5kZWZpbmVkLCBbc2xvd09yaWVudF0uY29uY2F0KENBQ0hFRCkpXG4gIGZvcih2YXIgaT0wOyBpPD1OVU1fRVhQQU5EOyArK2kpIHtcbiAgICBtb2R1bGUuZXhwb3J0c1tpXSA9IENBQ0hFRFtpXVxuICB9XG59XG5cbmdlbmVyYXRlT3JpZW50YXRpb25Qcm9jKCkiLCIvKipcbiAqIEJpdCB0d2lkZGxpbmcgaGFja3MgZm9yIEphdmFTY3JpcHQuXG4gKlxuICogQXV0aG9yOiBNaWtvbGEgTHlzZW5rb1xuICpcbiAqIFBvcnRlZCBmcm9tIFN0YW5mb3JkIGJpdCB0d2lkZGxpbmcgaGFjayBsaWJyYXJ5OlxuICogICAgaHR0cDovL2dyYXBoaWNzLnN0YW5mb3JkLmVkdS9+c2VhbmRlci9iaXRoYWNrcy5odG1sXG4gKi9cblxuXCJ1c2Ugc3RyaWN0XCI7IFwidXNlIHJlc3RyaWN0XCI7XG5cbi8vTnVtYmVyIG9mIGJpdHMgaW4gYW4gaW50ZWdlclxudmFyIElOVF9CSVRTID0gMzI7XG5cbi8vQ29uc3RhbnRzXG5leHBvcnRzLklOVF9CSVRTICA9IElOVF9CSVRTO1xuZXhwb3J0cy5JTlRfTUFYICAgPSAgMHg3ZmZmZmZmZjtcbmV4cG9ydHMuSU5UX01JTiAgID0gLTE8PChJTlRfQklUUy0xKTtcblxuLy9SZXR1cm5zIC0xLCAwLCArMSBkZXBlbmRpbmcgb24gc2lnbiBvZiB4XG5leHBvcnRzLnNpZ24gPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAodiA+IDApIC0gKHYgPCAwKTtcbn1cblxuLy9Db21wdXRlcyBhYnNvbHV0ZSB2YWx1ZSBvZiBpbnRlZ2VyXG5leHBvcnRzLmFicyA9IGZ1bmN0aW9uKHYpIHtcbiAgdmFyIG1hc2sgPSB2ID4+IChJTlRfQklUUy0xKTtcbiAgcmV0dXJuICh2IF4gbWFzaykgLSBtYXNrO1xufVxuXG4vL0NvbXB1dGVzIG1pbmltdW0gb2YgaW50ZWdlcnMgeCBhbmQgeVxuZXhwb3J0cy5taW4gPSBmdW5jdGlvbih4LCB5KSB7XG4gIHJldHVybiB5IF4gKCh4IF4geSkgJiAtKHggPCB5KSk7XG59XG5cbi8vQ29tcHV0ZXMgbWF4aW11bSBvZiBpbnRlZ2VycyB4IGFuZCB5XG5leHBvcnRzLm1heCA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgcmV0dXJuIHggXiAoKHggXiB5KSAmIC0oeCA8IHkpKTtcbn1cblxuLy9DaGVja3MgaWYgYSBudW1iZXIgaXMgYSBwb3dlciBvZiB0d29cbmV4cG9ydHMuaXNQb3cyID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gISh2ICYgKHYtMSkpICYmICghIXYpO1xufVxuXG4vL0NvbXB1dGVzIGxvZyBiYXNlIDIgb2YgdlxuZXhwb3J0cy5sb2cyID0gZnVuY3Rpb24odikge1xuICB2YXIgciwgc2hpZnQ7XG4gIHIgPSAgICAgKHYgPiAweEZGRkYpIDw8IDQ7IHYgPj4+PSByO1xuICBzaGlmdCA9ICh2ID4gMHhGRiAgKSA8PCAzOyB2ID4+Pj0gc2hpZnQ7IHIgfD0gc2hpZnQ7XG4gIHNoaWZ0ID0gKHYgPiAweEYgICApIDw8IDI7IHYgPj4+PSBzaGlmdDsgciB8PSBzaGlmdDtcbiAgc2hpZnQgPSAodiA+IDB4MyAgICkgPDwgMTsgdiA+Pj49IHNoaWZ0OyByIHw9IHNoaWZ0O1xuICByZXR1cm4gciB8ICh2ID4+IDEpO1xufVxuXG4vL0NvbXB1dGVzIGxvZyBiYXNlIDEwIG9mIHZcbmV4cG9ydHMubG9nMTAgPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAgKHYgPj0gMTAwMDAwMDAwMCkgPyA5IDogKHYgPj0gMTAwMDAwMDAwKSA/IDggOiAodiA+PSAxMDAwMDAwMCkgPyA3IDpcbiAgICAgICAgICAodiA+PSAxMDAwMDAwKSA/IDYgOiAodiA+PSAxMDAwMDApID8gNSA6ICh2ID49IDEwMDAwKSA/IDQgOlxuICAgICAgICAgICh2ID49IDEwMDApID8gMyA6ICh2ID49IDEwMCkgPyAyIDogKHYgPj0gMTApID8gMSA6IDA7XG59XG5cbi8vQ291bnRzIG51bWJlciBvZiBiaXRzXG5leHBvcnRzLnBvcENvdW50ID0gZnVuY3Rpb24odikge1xuICB2ID0gdiAtICgodiA+Pj4gMSkgJiAweDU1NTU1NTU1KTtcbiAgdiA9ICh2ICYgMHgzMzMzMzMzMykgKyAoKHYgPj4+IDIpICYgMHgzMzMzMzMzMyk7XG4gIHJldHVybiAoKHYgKyAodiA+Pj4gNCkgJiAweEYwRjBGMEYpICogMHgxMDEwMTAxKSA+Pj4gMjQ7XG59XG5cbi8vQ291bnRzIG51bWJlciBvZiB0cmFpbGluZyB6ZXJvc1xuZnVuY3Rpb24gY291bnRUcmFpbGluZ1plcm9zKHYpIHtcbiAgdmFyIGMgPSAzMjtcbiAgdiAmPSAtdjtcbiAgaWYgKHYpIGMtLTtcbiAgaWYgKHYgJiAweDAwMDBGRkZGKSBjIC09IDE2O1xuICBpZiAodiAmIDB4MDBGRjAwRkYpIGMgLT0gODtcbiAgaWYgKHYgJiAweDBGMEYwRjBGKSBjIC09IDQ7XG4gIGlmICh2ICYgMHgzMzMzMzMzMykgYyAtPSAyO1xuICBpZiAodiAmIDB4NTU1NTU1NTUpIGMgLT0gMTtcbiAgcmV0dXJuIGM7XG59XG5leHBvcnRzLmNvdW50VHJhaWxpbmdaZXJvcyA9IGNvdW50VHJhaWxpbmdaZXJvcztcblxuLy9Sb3VuZHMgdG8gbmV4dCBwb3dlciBvZiAyXG5leHBvcnRzLm5leHRQb3cyID0gZnVuY3Rpb24odikge1xuICB2ICs9IHYgPT09IDA7XG4gIC0tdjtcbiAgdiB8PSB2ID4+PiAxO1xuICB2IHw9IHYgPj4+IDI7XG4gIHYgfD0gdiA+Pj4gNDtcbiAgdiB8PSB2ID4+PiA4O1xuICB2IHw9IHYgPj4+IDE2O1xuICByZXR1cm4gdiArIDE7XG59XG5cbi8vUm91bmRzIGRvd24gdG8gcHJldmlvdXMgcG93ZXIgb2YgMlxuZXhwb3J0cy5wcmV2UG93MiA9IGZ1bmN0aW9uKHYpIHtcbiAgdiB8PSB2ID4+PiAxO1xuICB2IHw9IHYgPj4+IDI7XG4gIHYgfD0gdiA+Pj4gNDtcbiAgdiB8PSB2ID4+PiA4O1xuICB2IHw9IHYgPj4+IDE2O1xuICByZXR1cm4gdiAtICh2Pj4+MSk7XG59XG5cbi8vQ29tcHV0ZXMgcGFyaXR5IG9mIHdvcmRcbmV4cG9ydHMucGFyaXR5ID0gZnVuY3Rpb24odikge1xuICB2IF49IHYgPj4+IDE2O1xuICB2IF49IHYgPj4+IDg7XG4gIHYgXj0gdiA+Pj4gNDtcbiAgdiAmPSAweGY7XG4gIHJldHVybiAoMHg2OTk2ID4+PiB2KSAmIDE7XG59XG5cbnZhciBSRVZFUlNFX1RBQkxFID0gbmV3IEFycmF5KDI1Nik7XG5cbihmdW5jdGlvbih0YWIpIHtcbiAgZm9yKHZhciBpPTA7IGk8MjU2OyArK2kpIHtcbiAgICB2YXIgdiA9IGksIHIgPSBpLCBzID0gNztcbiAgICBmb3IgKHYgPj4+PSAxOyB2OyB2ID4+Pj0gMSkge1xuICAgICAgciA8PD0gMTtcbiAgICAgIHIgfD0gdiAmIDE7XG4gICAgICAtLXM7XG4gICAgfVxuICAgIHRhYltpXSA9IChyIDw8IHMpICYgMHhmZjtcbiAgfVxufSkoUkVWRVJTRV9UQUJMRSk7XG5cbi8vUmV2ZXJzZSBiaXRzIGluIGEgMzIgYml0IHdvcmRcbmV4cG9ydHMucmV2ZXJzZSA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICAoUkVWRVJTRV9UQUJMRVsgdiAgICAgICAgICYgMHhmZl0gPDwgMjQpIHxcbiAgICAgICAgICAoUkVWRVJTRV9UQUJMRVsodiA+Pj4gOCkgICYgMHhmZl0gPDwgMTYpIHxcbiAgICAgICAgICAoUkVWRVJTRV9UQUJMRVsodiA+Pj4gMTYpICYgMHhmZl0gPDwgOCkgIHxcbiAgICAgICAgICAgUkVWRVJTRV9UQUJMRVsodiA+Pj4gMjQpICYgMHhmZl07XG59XG5cbi8vSW50ZXJsZWF2ZSBiaXRzIG9mIDIgY29vcmRpbmF0ZXMgd2l0aCAxNiBiaXRzLiAgVXNlZnVsIGZvciBmYXN0IHF1YWR0cmVlIGNvZGVzXG5leHBvcnRzLmludGVybGVhdmUyID0gZnVuY3Rpb24oeCwgeSkge1xuICB4ICY9IDB4RkZGRjtcbiAgeCA9ICh4IHwgKHggPDwgOCkpICYgMHgwMEZGMDBGRjtcbiAgeCA9ICh4IHwgKHggPDwgNCkpICYgMHgwRjBGMEYwRjtcbiAgeCA9ICh4IHwgKHggPDwgMikpICYgMHgzMzMzMzMzMztcbiAgeCA9ICh4IHwgKHggPDwgMSkpICYgMHg1NTU1NTU1NTtcblxuICB5ICY9IDB4RkZGRjtcbiAgeSA9ICh5IHwgKHkgPDwgOCkpICYgMHgwMEZGMDBGRjtcbiAgeSA9ICh5IHwgKHkgPDwgNCkpICYgMHgwRjBGMEYwRjtcbiAgeSA9ICh5IHwgKHkgPDwgMikpICYgMHgzMzMzMzMzMztcbiAgeSA9ICh5IHwgKHkgPDwgMSkpICYgMHg1NTU1NTU1NTtcblxuICByZXR1cm4geCB8ICh5IDw8IDEpO1xufVxuXG4vL0V4dHJhY3RzIHRoZSBudGggaW50ZXJsZWF2ZWQgY29tcG9uZW50XG5leHBvcnRzLmRlaW50ZXJsZWF2ZTIgPSBmdW5jdGlvbih2LCBuKSB7XG4gIHYgPSAodiA+Pj4gbikgJiAweDU1NTU1NTU1O1xuICB2ID0gKHYgfCAodiA+Pj4gMSkpICAmIDB4MzMzMzMzMzM7XG4gIHYgPSAodiB8ICh2ID4+PiAyKSkgICYgMHgwRjBGMEYwRjtcbiAgdiA9ICh2IHwgKHYgPj4+IDQpKSAgJiAweDAwRkYwMEZGO1xuICB2ID0gKHYgfCAodiA+Pj4gMTYpKSAmIDB4MDAwRkZGRjtcbiAgcmV0dXJuICh2IDw8IDE2KSA+PiAxNjtcbn1cblxuXG4vL0ludGVybGVhdmUgYml0cyBvZiAzIGNvb3JkaW5hdGVzLCBlYWNoIHdpdGggMTAgYml0cy4gIFVzZWZ1bCBmb3IgZmFzdCBvY3RyZWUgY29kZXNcbmV4cG9ydHMuaW50ZXJsZWF2ZTMgPSBmdW5jdGlvbih4LCB5LCB6KSB7XG4gIHggJj0gMHgzRkY7XG4gIHggID0gKHggfCAoeDw8MTYpKSAmIDQyNzgxOTAzMzU7XG4gIHggID0gKHggfCAoeDw8OCkpICAmIDI1MTcxOTY5NTtcbiAgeCAgPSAoeCB8ICh4PDw0KSkgICYgMzI3MjM1NjAzNTtcbiAgeCAgPSAoeCB8ICh4PDwyKSkgICYgMTIyNzEzMzUxMztcblxuICB5ICY9IDB4M0ZGO1xuICB5ICA9ICh5IHwgKHk8PDE2KSkgJiA0Mjc4MTkwMzM1O1xuICB5ICA9ICh5IHwgKHk8PDgpKSAgJiAyNTE3MTk2OTU7XG4gIHkgID0gKHkgfCAoeTw8NCkpICAmIDMyNzIzNTYwMzU7XG4gIHkgID0gKHkgfCAoeTw8MikpICAmIDEyMjcxMzM1MTM7XG4gIHggfD0gKHkgPDwgMSk7XG4gIFxuICB6ICY9IDB4M0ZGO1xuICB6ICA9ICh6IHwgKHo8PDE2KSkgJiA0Mjc4MTkwMzM1O1xuICB6ICA9ICh6IHwgKHo8PDgpKSAgJiAyNTE3MTk2OTU7XG4gIHogID0gKHogfCAoejw8NCkpICAmIDMyNzIzNTYwMzU7XG4gIHogID0gKHogfCAoejw8MikpICAmIDEyMjcxMzM1MTM7XG4gIFxuICByZXR1cm4geCB8ICh6IDw8IDIpO1xufVxuXG4vL0V4dHJhY3RzIG50aCBpbnRlcmxlYXZlZCBjb21wb25lbnQgb2YgYSAzLXR1cGxlXG5leHBvcnRzLmRlaW50ZXJsZWF2ZTMgPSBmdW5jdGlvbih2LCBuKSB7XG4gIHYgPSAodiA+Pj4gbikgICAgICAgJiAxMjI3MTMzNTEzO1xuICB2ID0gKHYgfCAodj4+PjIpKSAgICYgMzI3MjM1NjAzNTtcbiAgdiA9ICh2IHwgKHY+Pj40KSkgICAmIDI1MTcxOTY5NTtcbiAgdiA9ICh2IHwgKHY+Pj44KSkgICAmIDQyNzgxOTAzMzU7XG4gIHYgPSAodiB8ICh2Pj4+MTYpKSAgJiAweDNGRjtcbiAgcmV0dXJuICh2PDwyMik+PjIyO1xufVxuXG4vL0NvbXB1dGVzIG5leHQgY29tYmluYXRpb24gaW4gY29sZXhpY29ncmFwaGljIG9yZGVyICh0aGlzIGlzIG1pc3Rha2VubHkgY2FsbGVkIG5leHRQZXJtdXRhdGlvbiBvbiB0aGUgYml0IHR3aWRkbGluZyBoYWNrcyBwYWdlKVxuZXhwb3J0cy5uZXh0Q29tYmluYXRpb24gPSBmdW5jdGlvbih2KSB7XG4gIHZhciB0ID0gdiB8ICh2IC0gMSk7XG4gIHJldHVybiAodCArIDEpIHwgKCgofnQgJiAtfnQpIC0gMSkgPj4+IChjb3VudFRyYWlsaW5nWmVyb3ModikgKyAxKSk7XG59XG5cbiIsIlwidXNlIHN0cmljdFwiOyBcInVzZSByZXN0cmljdFwiO1xuXG5tb2R1bGUuZXhwb3J0cyA9IFVuaW9uRmluZDtcblxuZnVuY3Rpb24gVW5pb25GaW5kKGNvdW50KSB7XG4gIHRoaXMucm9vdHMgPSBuZXcgQXJyYXkoY291bnQpO1xuICB0aGlzLnJhbmtzID0gbmV3IEFycmF5KGNvdW50KTtcbiAgXG4gIGZvcih2YXIgaT0wOyBpPGNvdW50OyArK2kpIHtcbiAgICB0aGlzLnJvb3RzW2ldID0gaTtcbiAgICB0aGlzLnJhbmtzW2ldID0gMDtcbiAgfVxufVxuXG52YXIgcHJvdG8gPSBVbmlvbkZpbmQucHJvdG90eXBlXG5cbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShwcm90bywgXCJsZW5ndGhcIiwge1xuICBcImdldFwiOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5yb290cy5sZW5ndGhcbiAgfVxufSlcblxucHJvdG8ubWFrZVNldCA9IGZ1bmN0aW9uKCkge1xuICB2YXIgbiA9IHRoaXMucm9vdHMubGVuZ3RoO1xuICB0aGlzLnJvb3RzLnB1c2gobik7XG4gIHRoaXMucmFua3MucHVzaCgwKTtcbiAgcmV0dXJuIG47XG59XG5cbnByb3RvLmZpbmQgPSBmdW5jdGlvbih4KSB7XG4gIHZhciByb290cyA9IHRoaXMucm9vdHM7XG4gIHdoaWxlKHJvb3RzW3hdICE9PSB4KSB7XG4gICAgdmFyIHkgPSByb290c1t4XTtcbiAgICByb290c1t4XSA9IHJvb3RzW3ldO1xuICAgIHggPSB5O1xuICB9XG4gIHJldHVybiB4O1xufVxuXG5wcm90by5saW5rID0gZnVuY3Rpb24oeCwgeSkge1xuICB2YXIgeHIgPSB0aGlzLmZpbmQoeClcbiAgICAsIHlyID0gdGhpcy5maW5kKHkpO1xuICBpZih4ciA9PT0geXIpIHtcbiAgICByZXR1cm47XG4gIH1cbiAgdmFyIHJhbmtzID0gdGhpcy5yYW5rc1xuICAgICwgcm9vdHMgPSB0aGlzLnJvb3RzXG4gICAgLCB4ZCAgICA9IHJhbmtzW3hyXVxuICAgICwgeWQgICAgPSByYW5rc1t5cl07XG4gIGlmKHhkIDwgeWQpIHtcbiAgICByb290c1t4cl0gPSB5cjtcbiAgfSBlbHNlIGlmKHlkIDwgeGQpIHtcbiAgICByb290c1t5cl0gPSB4cjtcbiAgfSBlbHNlIHtcbiAgICByb290c1t5cl0gPSB4cjtcbiAgICArK3JhbmtzW3hyXTtcbiAgfVxufSIsIlwidXNlIHN0cmljdFwiOyBcInVzZSByZXN0cmljdFwiO1xuXG52YXIgYml0cyAgICAgID0gcmVxdWlyZShcImJpdC10d2lkZGxlXCIpXG4gICwgVW5pb25GaW5kID0gcmVxdWlyZShcInVuaW9uLWZpbmRcIilcblxuLy9SZXR1cm5zIHRoZSBkaW1lbnNpb24gb2YgYSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGRpbWVuc2lvbihjZWxscykge1xuICB2YXIgZCA9IDBcbiAgICAsIG1heCA9IE1hdGgubWF4XG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIGQgPSBtYXgoZCwgY2VsbHNbaV0ubGVuZ3RoKVxuICB9XG4gIHJldHVybiBkLTFcbn1cbmV4cG9ydHMuZGltZW5zaW9uID0gZGltZW5zaW9uXG5cbi8vQ291bnRzIHRoZSBudW1iZXIgb2YgdmVydGljZXMgaW4gZmFjZXNcbmZ1bmN0aW9uIGNvdW50VmVydGljZXMoY2VsbHMpIHtcbiAgdmFyIHZjID0gLTFcbiAgICAsIG1heCA9IE1hdGgubWF4XG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MCwgamw9Yy5sZW5ndGg7IGo8amw7ICsraikge1xuICAgICAgdmMgPSBtYXgodmMsIGNbal0pXG4gICAgfVxuICB9XG4gIHJldHVybiB2YysxXG59XG5leHBvcnRzLmNvdW50VmVydGljZXMgPSBjb3VudFZlcnRpY2VzXG5cbi8vUmV0dXJucyBhIGRlZXAgY29weSBvZiBjZWxsc1xuZnVuY3Rpb24gY2xvbmVDZWxscyhjZWxscykge1xuICB2YXIgbmNlbGxzID0gbmV3IEFycmF5KGNlbGxzLmxlbmd0aClcbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgbmNlbGxzW2ldID0gY2VsbHNbaV0uc2xpY2UoMClcbiAgfVxuICByZXR1cm4gbmNlbGxzXG59XG5leHBvcnRzLmNsb25lQ2VsbHMgPSBjbG9uZUNlbGxzXG5cbi8vUmFua3MgYSBwYWlyIG9mIGNlbGxzIHVwIHRvIHBlcm11dGF0aW9uXG5mdW5jdGlvbiBjb21wYXJlQ2VsbHMoYSwgYikge1xuICB2YXIgbiA9IGEubGVuZ3RoXG4gICAgLCB0ID0gYS5sZW5ndGggLSBiLmxlbmd0aFxuICAgICwgbWluID0gTWF0aC5taW5cbiAgaWYodCkge1xuICAgIHJldHVybiB0XG4gIH1cbiAgc3dpdGNoKG4pIHtcbiAgICBjYXNlIDA6XG4gICAgICByZXR1cm4gMDtcbiAgICBjYXNlIDE6XG4gICAgICByZXR1cm4gYVswXSAtIGJbMF07XG4gICAgY2FzZSAyOlxuICAgICAgdmFyIGQgPSBhWzBdK2FbMV0tYlswXS1iWzFdXG4gICAgICBpZihkKSB7XG4gICAgICAgIHJldHVybiBkXG4gICAgICB9XG4gICAgICByZXR1cm4gbWluKGFbMF0sYVsxXSkgLSBtaW4oYlswXSxiWzFdKVxuICAgIGNhc2UgMzpcbiAgICAgIHZhciBsMSA9IGFbMF0rYVsxXVxuICAgICAgICAsIG0xID0gYlswXStiWzFdXG4gICAgICBkID0gbDErYVsyXSAtIChtMStiWzJdKVxuICAgICAgaWYoZCkge1xuICAgICAgICByZXR1cm4gZFxuICAgICAgfVxuICAgICAgdmFyIGwwID0gbWluKGFbMF0sIGFbMV0pXG4gICAgICAgICwgbTAgPSBtaW4oYlswXSwgYlsxXSlcbiAgICAgICAgLCBkICA9IG1pbihsMCwgYVsyXSkgLSBtaW4obTAsIGJbMl0pXG4gICAgICBpZihkKSB7XG4gICAgICAgIHJldHVybiBkXG4gICAgICB9XG4gICAgICByZXR1cm4gbWluKGwwK2FbMl0sIGwxKSAtIG1pbihtMCtiWzJdLCBtMSlcbiAgICBcbiAgICAvL1RPRE86IE1heWJlIG9wdGltaXplIG49NCBhcyB3ZWxsP1xuICAgIFxuICAgIGRlZmF1bHQ6XG4gICAgICB2YXIgYXMgPSBhLnNsaWNlKDApXG4gICAgICBhcy5zb3J0KClcbiAgICAgIHZhciBicyA9IGIuc2xpY2UoMClcbiAgICAgIGJzLnNvcnQoKVxuICAgICAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgICAgIHQgPSBhc1tpXSAtIGJzW2ldXG4gICAgICAgIGlmKHQpIHtcbiAgICAgICAgICByZXR1cm4gdFxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXR1cm4gMFxuICB9XG59XG5leHBvcnRzLmNvbXBhcmVDZWxscyA9IGNvbXBhcmVDZWxsc1xuXG5mdW5jdGlvbiBjb21wYXJlWmlwcGVkKGEsIGIpIHtcbiAgcmV0dXJuIGNvbXBhcmVDZWxscyhhWzBdLCBiWzBdKVxufVxuXG4vL1B1dHMgYSBjZWxsIGNvbXBsZXggaW50byBub3JtYWwgb3JkZXIgZm9yIHRoZSBwdXJwb3NlcyBvZiBmaW5kQ2VsbCBxdWVyaWVzXG5mdW5jdGlvbiBub3JtYWxpemUoY2VsbHMsIGF0dHIpIHtcbiAgaWYoYXR0cikge1xuICAgIHZhciBsZW4gPSBjZWxscy5sZW5ndGhcbiAgICB2YXIgemlwcGVkID0gbmV3IEFycmF5KGxlbilcbiAgICBmb3IodmFyIGk9MDsgaTxsZW47ICsraSkge1xuICAgICAgemlwcGVkW2ldID0gW2NlbGxzW2ldLCBhdHRyW2ldXVxuICAgIH1cbiAgICB6aXBwZWQuc29ydChjb21wYXJlWmlwcGVkKVxuICAgIGZvcih2YXIgaT0wOyBpPGxlbjsgKytpKSB7XG4gICAgICBjZWxsc1tpXSA9IHppcHBlZFtpXVswXVxuICAgICAgYXR0cltpXSA9IHppcHBlZFtpXVsxXVxuICAgIH1cbiAgICByZXR1cm4gY2VsbHNcbiAgfSBlbHNlIHtcbiAgICBjZWxscy5zb3J0KGNvbXBhcmVDZWxscylcbiAgICByZXR1cm4gY2VsbHNcbiAgfVxufVxuZXhwb3J0cy5ub3JtYWxpemUgPSBub3JtYWxpemVcblxuLy9SZW1vdmVzIGFsbCBkdXBsaWNhdGUgY2VsbHMgaW4gdGhlIGNvbXBsZXhcbmZ1bmN0aW9uIHVuaXF1ZShjZWxscykge1xuICBpZihjZWxscy5sZW5ndGggPT09IDApIHtcbiAgICByZXR1cm4gW11cbiAgfVxuICB2YXIgcHRyID0gMVxuICAgICwgbGVuID0gY2VsbHMubGVuZ3RoXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpKSB7XG4gICAgdmFyIGEgPSBjZWxsc1tpXVxuICAgIGlmKGNvbXBhcmVDZWxscyhhLCBjZWxsc1tpLTFdKSkge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBjZWxsc1twdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGNlbGxzLmxlbmd0aCA9IHB0clxuICByZXR1cm4gY2VsbHNcbn1cbmV4cG9ydHMudW5pcXVlID0gdW5pcXVlO1xuXG4vL0ZpbmRzIGEgY2VsbCBpbiBhIG5vcm1hbGl6ZWQgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBmaW5kQ2VsbChjZWxscywgYykge1xuICB2YXIgbG8gPSAwXG4gICAgLCBoaSA9IGNlbGxzLmxlbmd0aC0xXG4gICAgLCByICA9IC0xXG4gIHdoaWxlIChsbyA8PSBoaSkge1xuICAgIHZhciBtaWQgPSAobG8gKyBoaSkgPj4gMVxuICAgICAgLCBzICAgPSBjb21wYXJlQ2VsbHMoY2VsbHNbbWlkXSwgYylcbiAgICBpZihzIDw9IDApIHtcbiAgICAgIGlmKHMgPT09IDApIHtcbiAgICAgICAgciA9IG1pZFxuICAgICAgfVxuICAgICAgbG8gPSBtaWQgKyAxXG4gICAgfSBlbHNlIGlmKHMgPiAwKSB7XG4gICAgICBoaSA9IG1pZCAtIDFcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJcbn1cbmV4cG9ydHMuZmluZENlbGwgPSBmaW5kQ2VsbDtcblxuLy9CdWlsZHMgYW4gaW5kZXggZm9yIGFuIG4tY2VsbC4gIFRoaXMgaXMgbW9yZSBnZW5lcmFsIHRoYW4gZHVhbCwgYnV0IGxlc3MgZWZmaWNpZW50XG5mdW5jdGlvbiBpbmNpZGVuY2UoZnJvbV9jZWxscywgdG9fY2VsbHMpIHtcbiAgdmFyIGluZGV4ID0gbmV3IEFycmF5KGZyb21fY2VsbHMubGVuZ3RoKVxuICBmb3IodmFyIGk9MCwgaWw9aW5kZXgubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICBpbmRleFtpXSA9IFtdXG4gIH1cbiAgdmFyIGIgPSBbXVxuICBmb3IodmFyIGk9MCwgbj10b19jZWxscy5sZW5ndGg7IGk8bjsgKytpKSB7XG4gICAgdmFyIGMgPSB0b19jZWxsc1tpXVxuICAgIHZhciBjbCA9IGMubGVuZ3RoXG4gICAgZm9yKHZhciBrPTEsIGtuPSgxPDxjbCk7IGs8a247ICsraykge1xuICAgICAgYi5sZW5ndGggPSBiaXRzLnBvcENvdW50KGspXG4gICAgICB2YXIgbCA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPGNsOyArK2opIHtcbiAgICAgICAgaWYoayAmICgxPDxqKSkge1xuICAgICAgICAgIGJbbCsrXSA9IGNbal1cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgdmFyIGlkeD1maW5kQ2VsbChmcm9tX2NlbGxzLCBiKVxuICAgICAgaWYoaWR4IDwgMCkge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgd2hpbGUodHJ1ZSkge1xuICAgICAgICBpbmRleFtpZHgrK10ucHVzaChpKVxuICAgICAgICBpZihpZHggPj0gZnJvbV9jZWxscy5sZW5ndGggfHwgY29tcGFyZUNlbGxzKGZyb21fY2VsbHNbaWR4XSwgYikgIT09IDApIHtcbiAgICAgICAgICBicmVha1xuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICB9XG4gIHJldHVybiBpbmRleFxufVxuZXhwb3J0cy5pbmNpZGVuY2UgPSBpbmNpZGVuY2VcblxuLy9Db21wdXRlcyB0aGUgZHVhbCBvZiB0aGUgbWVzaC4gIFRoaXMgaXMgYmFzaWNhbGx5IGFuIG9wdGltaXplZCB2ZXJzaW9uIG9mIGJ1aWxkSW5kZXggZm9yIHRoZSBzaXR1YXRpb24gd2hlcmUgZnJvbV9jZWxscyBpcyBqdXN0IHRoZSBsaXN0IG9mIHZlcnRpY2VzXG5mdW5jdGlvbiBkdWFsKGNlbGxzLCB2ZXJ0ZXhfY291bnQpIHtcbiAgaWYoIXZlcnRleF9jb3VudCkge1xuICAgIHJldHVybiBpbmNpZGVuY2UodW5pcXVlKHNrZWxldG9uKGNlbGxzLCAwKSksIGNlbGxzLCAwKVxuICB9XG4gIHZhciByZXMgPSBuZXcgQXJyYXkodmVydGV4X2NvdW50KVxuICBmb3IodmFyIGk9MDsgaTx2ZXJ0ZXhfY291bnQ7ICsraSkge1xuICAgIHJlc1tpXSA9IFtdXG4gIH1cbiAgZm9yKHZhciBpPTAsIGxlbj1jZWxscy5sZW5ndGg7IGk8bGVuOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTAsIGNsPWMubGVuZ3RoOyBqPGNsOyArK2opIHtcbiAgICAgIHJlc1tjW2pdXS5wdXNoKGkpXG4gICAgfVxuICB9XG4gIHJldHVybiByZXNcbn1cbmV4cG9ydHMuZHVhbCA9IGR1YWxcblxuLy9FbnVtZXJhdGVzIGFsbCBjZWxscyBpbiB0aGUgY29tcGxleFxuZnVuY3Rpb24gZXhwbG9kZShjZWxscykge1xuICB2YXIgcmVzdWx0ID0gW11cbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgICAgLCBjbCA9IGMubGVuZ3RofDBcbiAgICBmb3IodmFyIGo9MSwgamw9KDE8PGNsKTsgajxqbDsgKytqKSB7XG4gICAgICB2YXIgYiA9IFtdXG4gICAgICBmb3IodmFyIGs9MDsgazxjbDsgKytrKSB7XG4gICAgICAgIGlmKChqID4+PiBrKSAmIDEpIHtcbiAgICAgICAgICBiLnB1c2goY1trXSlcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmVzdWx0LnB1c2goYilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIG5vcm1hbGl6ZShyZXN1bHQpXG59XG5leHBvcnRzLmV4cGxvZGUgPSBleHBsb2RlXG5cbi8vRW51bWVyYXRlcyBhbGwgb2YgdGhlIG4tY2VsbHMgb2YgYSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIHNrZWxldG9uKGNlbGxzLCBuKSB7XG4gIGlmKG4gPCAwKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cbiAgdmFyIHJlc3VsdCA9IFtdXG4gICAgLCBrMCAgICAgPSAoMTw8KG4rMSkpLTFcbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBrPWswOyBrPCgxPDxjLmxlbmd0aCk7IGs9Yml0cy5uZXh0Q29tYmluYXRpb24oaykpIHtcbiAgICAgIHZhciBiID0gbmV3IEFycmF5KG4rMSlcbiAgICAgICAgLCBsID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8Yy5sZW5ndGg7ICsraikge1xuICAgICAgICBpZihrICYgKDE8PGopKSB7XG4gICAgICAgICAgYltsKytdID0gY1tqXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXN1bHQucHVzaChiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gbm9ybWFsaXplKHJlc3VsdClcbn1cbmV4cG9ydHMuc2tlbGV0b24gPSBza2VsZXRvbjtcblxuLy9Db21wdXRlcyB0aGUgYm91bmRhcnkgb2YgYWxsIGNlbGxzLCBkb2VzIG5vdCByZW1vdmUgZHVwbGljYXRlc1xuZnVuY3Rpb24gYm91bmRhcnkoY2VsbHMpIHtcbiAgdmFyIHJlcyA9IFtdXG4gIGZvcih2YXIgaT0wLGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wLGNsPWMubGVuZ3RoOyBqPGNsOyArK2opIHtcbiAgICAgIHZhciBiID0gbmV3IEFycmF5KGMubGVuZ3RoLTEpXG4gICAgICBmb3IodmFyIGs9MCwgbD0wOyBrPGNsOyArK2spIHtcbiAgICAgICAgaWYoayAhPT0gaikge1xuICAgICAgICAgIGJbbCsrXSA9IGNba11cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmVzLnB1c2goYilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIG5vcm1hbGl6ZShyZXMpXG59XG5leHBvcnRzLmJvdW5kYXJ5ID0gYm91bmRhcnk7XG5cbi8vQ29tcHV0ZXMgY29ubmVjdGVkIGNvbXBvbmVudHMgZm9yIGEgZGVuc2UgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBjb25uZWN0ZWRDb21wb25lbnRzX2RlbnNlKGNlbGxzLCB2ZXJ0ZXhfY291bnQpIHtcbiAgdmFyIGxhYmVscyA9IG5ldyBVbmlvbkZpbmQodmVydGV4X2NvdW50KVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MDsgajxjLmxlbmd0aDsgKytqKSB7XG4gICAgICBmb3IodmFyIGs9aisxOyBrPGMubGVuZ3RoOyArK2spIHtcbiAgICAgICAgbGFiZWxzLmxpbmsoY1tqXSwgY1trXSlcbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgdmFyIGNvbXBvbmVudHMgPSBbXVxuICAgICwgY29tcG9uZW50X2xhYmVscyA9IGxhYmVscy5yYW5rc1xuICBmb3IodmFyIGk9MDsgaTxjb21wb25lbnRfbGFiZWxzLmxlbmd0aDsgKytpKSB7XG4gICAgY29tcG9uZW50X2xhYmVsc1tpXSA9IC0xXG4gIH1cbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgbCA9IGxhYmVscy5maW5kKGNlbGxzW2ldWzBdKVxuICAgIGlmKGNvbXBvbmVudF9sYWJlbHNbbF0gPCAwKSB7XG4gICAgICBjb21wb25lbnRfbGFiZWxzW2xdID0gY29tcG9uZW50cy5sZW5ndGhcbiAgICAgIGNvbXBvbmVudHMucHVzaChbY2VsbHNbaV0uc2xpY2UoMCldKVxuICAgIH0gZWxzZSB7XG4gICAgICBjb21wb25lbnRzW2NvbXBvbmVudF9sYWJlbHNbbF1dLnB1c2goY2VsbHNbaV0uc2xpY2UoMCkpXG4gICAgfVxuICB9XG4gIHJldHVybiBjb21wb25lbnRzXG59XG5cbi8vQ29tcHV0ZXMgY29ubmVjdGVkIGNvbXBvbmVudHMgZm9yIGEgc3BhcnNlIGdyYXBoXG5mdW5jdGlvbiBjb25uZWN0ZWRDb21wb25lbnRzX3NwYXJzZShjZWxscykge1xuICB2YXIgdmVydGljZXMgID0gdW5pcXVlKG5vcm1hbGl6ZShza2VsZXRvbihjZWxscywgMCkpKVxuICAgICwgbGFiZWxzICAgID0gbmV3IFVuaW9uRmluZCh2ZXJ0aWNlcy5sZW5ndGgpXG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wOyBqPGMubGVuZ3RoOyArK2opIHtcbiAgICAgIHZhciB2aiA9IGZpbmRDZWxsKHZlcnRpY2VzLCBbY1tqXV0pXG4gICAgICBmb3IodmFyIGs9aisxOyBrPGMubGVuZ3RoOyArK2spIHtcbiAgICAgICAgbGFiZWxzLmxpbmsodmosIGZpbmRDZWxsKHZlcnRpY2VzLCBbY1trXV0pKVxuICAgICAgfVxuICAgIH1cbiAgfVxuICB2YXIgY29tcG9uZW50cyAgICAgICAgPSBbXVxuICAgICwgY29tcG9uZW50X2xhYmVscyAgPSBsYWJlbHMucmFua3NcbiAgZm9yKHZhciBpPTA7IGk8Y29tcG9uZW50X2xhYmVscy5sZW5ndGg7ICsraSkge1xuICAgIGNvbXBvbmVudF9sYWJlbHNbaV0gPSAtMVxuICB9XG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGwgPSBsYWJlbHMuZmluZChmaW5kQ2VsbCh2ZXJ0aWNlcywgW2NlbGxzW2ldWzBdXSkpO1xuICAgIGlmKGNvbXBvbmVudF9sYWJlbHNbbF0gPCAwKSB7XG4gICAgICBjb21wb25lbnRfbGFiZWxzW2xdID0gY29tcG9uZW50cy5sZW5ndGhcbiAgICAgIGNvbXBvbmVudHMucHVzaChbY2VsbHNbaV0uc2xpY2UoMCldKVxuICAgIH0gZWxzZSB7XG4gICAgICBjb21wb25lbnRzW2NvbXBvbmVudF9sYWJlbHNbbF1dLnB1c2goY2VsbHNbaV0uc2xpY2UoMCkpXG4gICAgfVxuICB9XG4gIHJldHVybiBjb21wb25lbnRzXG59XG5cbi8vQ29tcHV0ZXMgY29ubmVjdGVkIGNvbXBvbmVudHMgZm9yIGEgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBjb25uZWN0ZWRDb21wb25lbnRzKGNlbGxzLCB2ZXJ0ZXhfY291bnQpIHtcbiAgaWYodmVydGV4X2NvdW50KSB7XG4gICAgcmV0dXJuIGNvbm5lY3RlZENvbXBvbmVudHNfZGVuc2UoY2VsbHMsIHZlcnRleF9jb3VudClcbiAgfVxuICByZXR1cm4gY29ubmVjdGVkQ29tcG9uZW50c19zcGFyc2UoY2VsbHMpXG59XG5leHBvcnRzLmNvbm5lY3RlZENvbXBvbmVudHMgPSBjb25uZWN0ZWRDb21wb25lbnRzXG4iLCJcInVzZSBzdHJpY3RcIlxuXG5mdW5jdGlvbiB1bmlxdWVfcHJlZChsaXN0LCBjb21wYXJlKSB7XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBsaXN0Lmxlbmd0aFxuICAgICwgYT1saXN0WzBdLCBiPWxpc3RbMF1cbiAgZm9yKHZhciBpPTE7IGk8bGVuOyArK2kpIHtcbiAgICBiID0gYVxuICAgIGEgPSBsaXN0W2ldXG4gICAgaWYoY29tcGFyZShhLCBiKSkge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBsaXN0W3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgbGlzdC5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGxpc3Rcbn1cblxuZnVuY3Rpb24gdW5pcXVlX2VxKGxpc3QpIHtcbiAgdmFyIHB0ciA9IDFcbiAgICAsIGxlbiA9IGxpc3QubGVuZ3RoXG4gICAgLCBhPWxpc3RbMF0sIGIgPSBsaXN0WzBdXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpLCBiPWEpIHtcbiAgICBiID0gYVxuICAgIGEgPSBsaXN0W2ldXG4gICAgaWYoYSAhPT0gYikge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBsaXN0W3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgbGlzdC5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGxpc3Rcbn1cblxuZnVuY3Rpb24gdW5pcXVlKGxpc3QsIGNvbXBhcmUsIHNvcnRlZCkge1xuICBpZihsaXN0Lmxlbmd0aCA9PT0gMCkge1xuICAgIHJldHVybiBsaXN0XG4gIH1cbiAgaWYoY29tcGFyZSkge1xuICAgIGlmKCFzb3J0ZWQpIHtcbiAgICAgIGxpc3Quc29ydChjb21wYXJlKVxuICAgIH1cbiAgICByZXR1cm4gdW5pcXVlX3ByZWQobGlzdCwgY29tcGFyZSlcbiAgfVxuICBpZighc29ydGVkKSB7XG4gICAgbGlzdC5zb3J0KClcbiAgfVxuICByZXR1cm4gdW5pcXVlX2VxKGxpc3QpXG59XG5cbm1vZHVsZS5leHBvcnRzID0gdW5pcXVlXG4iLCJcInVzZSBzdHJpY3RcIlxuXG52YXIgY2ggPSByZXF1aXJlKFwiaW5jcmVtZW50YWwtY29udmV4LWh1bGxcIilcbnZhciB1bmlxID0gcmVxdWlyZShcInVuaXFcIilcblxubW9kdWxlLmV4cG9ydHMgPSB0cmlhbmd1bGF0ZVxuXG5mdW5jdGlvbiBMaWZ0ZWRQb2ludChwLCBpKSB7XG4gIHRoaXMucG9pbnQgPSBwXG4gIHRoaXMuaW5kZXggPSBpXG59XG5cbmZ1bmN0aW9uIGNvbXBhcmVMaWZ0ZWQoYSwgYikge1xuICB2YXIgYXAgPSBhLnBvaW50XG4gIHZhciBicCA9IGIucG9pbnRcbiAgdmFyIGQgPSBhcC5sZW5ndGhcbiAgZm9yKHZhciBpPTA7IGk8ZDsgKytpKSB7XG4gICAgdmFyIHMgPSBicFtpXSAtIGFwW2ldXG4gICAgaWYocykge1xuICAgICAgcmV0dXJuIHNcbiAgICB9XG4gIH1cbiAgcmV0dXJuIDBcbn1cblxuZnVuY3Rpb24gdHJpYW5ndWxhdGUxRChuLCBwb2ludHMsIGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgaWYobiA9PT0gMSkge1xuICAgIGlmKGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgICAgIHJldHVybiBbIFstMSwgMF0gXVxuICAgIH0gZWxzZSB7XG4gICAgICByZXR1cm4gW11cbiAgICB9XG4gIH1cbiAgdmFyIGxpZnRlZCA9IHBvaW50cy5tYXAoZnVuY3Rpb24ocCwgaSkge1xuICAgIHJldHVybiBbIHBbMF0sIGkgXVxuICB9KVxuICBsaWZ0ZWQuc29ydChmdW5jdGlvbihhLGIpIHtcbiAgICByZXR1cm4gYVswXSAtIGJbMF1cbiAgfSlcbiAgdmFyIGNlbGxzID0gbmV3IEFycmF5KG4gLSAxKVxuICBmb3IodmFyIGk9MTsgaTxuOyArK2kpIHtcbiAgICB2YXIgYSA9IGxpZnRlZFtpLTFdXG4gICAgdmFyIGIgPSBsaWZ0ZWRbaV1cbiAgICBjZWxsc1tpLTFdID0gWyBhWzFdLCBiWzFdIF1cbiAgfVxuICBpZihpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gICAgY2VsbHMucHVzaChcbiAgICAgIFsgLTEsIGNlbGxzWzBdWzFdLCBdLFxuICAgICAgWyBjZWxsc1tuLTFdWzFdLCAtMSBdKVxuICB9XG4gIHJldHVybiBjZWxsc1xufVxuXG5mdW5jdGlvbiB0cmlhbmd1bGF0ZShwb2ludHMsIGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgdmFyIG4gPSBwb2ludHMubGVuZ3RoXG4gIGlmKG4gPT09IDApIHtcbiAgICByZXR1cm4gW11cbiAgfVxuICBcbiAgdmFyIGQgPSBwb2ludHNbMF0ubGVuZ3RoXG4gIGlmKGQgPCAxKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cblxuICAvL1NwZWNpYWwgY2FzZTogIEZvciAxRCB3ZSBjYW4ganVzdCBzb3J0IHRoZSBwb2ludHNcbiAgaWYoZCA9PT0gMSkge1xuICAgIHJldHVybiB0cmlhbmd1bGF0ZTFEKG4sIHBvaW50cywgaW5jbHVkZVBvaW50QXRJbmZpbml0eSlcbiAgfVxuICBcbiAgLy9MaWZ0IHBvaW50cywgc29ydFxuICB2YXIgbGlmdGVkID0gbmV3IEFycmF5KG4pXG4gIHZhciB1cHBlciA9IDEuMFxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICB2YXIgcCA9IHBvaW50c1tpXVxuICAgIHZhciB4ID0gbmV3IEFycmF5KGQrMSlcbiAgICB2YXIgbCA9IDAuMFxuICAgIGZvcih2YXIgaj0wOyBqPGQ7ICsraikge1xuICAgICAgdmFyIHYgPSBwW2pdXG4gICAgICB4W2pdID0gdlxuICAgICAgbCArPSB2ICogdlxuICAgIH1cbiAgICB4W2RdID0gbFxuICAgIGxpZnRlZFtpXSA9IG5ldyBMaWZ0ZWRQb2ludCh4LCBpKVxuICAgIHVwcGVyID0gTWF0aC5tYXgobCwgdXBwZXIpXG4gIH1cbiAgdW5pcShsaWZ0ZWQsIGNvbXBhcmVMaWZ0ZWQpXG4gIFxuICAvL0RvdWJsZSBwb2ludHNcbiAgbiA9IGxpZnRlZC5sZW5ndGhcblxuICAvL0NyZWF0ZSBuZXcgbGlzdCBvZiBwb2ludHNcbiAgdmFyIGRwb2ludHMgPSBuZXcgQXJyYXkobiArIGQgKyAxKVxuICB2YXIgZGluZGV4ID0gbmV3IEFycmF5KG4gKyBkICsgMSlcblxuICAvL0FkZCBzdGVpbmVyIHBvaW50cyBhdCB0b3BcbiAgdmFyIHUgPSAoZCsxKSAqIChkKzEpICogdXBwZXJcbiAgdmFyIHkgPSBuZXcgQXJyYXkoZCsxKVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgeVtpXSA9IDAuMFxuICB9XG4gIHlbZF0gPSB1XG5cbiAgZHBvaW50c1swXSA9IHkuc2xpY2UoKVxuICBkaW5kZXhbMF0gPSAtMVxuXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgeCA9IHkuc2xpY2UoKVxuICAgIHhbaV0gPSAxXG4gICAgZHBvaW50c1tpKzFdID0geFxuICAgIGRpbmRleFtpKzFdID0gLTFcbiAgfVxuXG4gIC8vQ29weSByZXN0IG9mIHRoZSBwb2ludHMgb3ZlclxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICB2YXIgaCA9IGxpZnRlZFtpXVxuICAgIGRwb2ludHNbaSArIGQgKyAxXSA9IGgucG9pbnRcbiAgICBkaW5kZXhbaSArIGQgKyAxXSA9ICBoLmluZGV4XG4gIH1cblxuICAvL0NvbnN0cnVjdCBjb252ZXggaHVsbFxuICB2YXIgaHVsbCA9IGNoKGRwb2ludHMsIGZhbHNlKVxuICBpZihpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gICAgaHVsbCA9IGh1bGwuZmlsdGVyKGZ1bmN0aW9uKGNlbGwpIHtcbiAgICAgIHZhciBjb3VudCA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgdmFyIHYgPSBkaW5kZXhbY2VsbFtqXV1cbiAgICAgICAgaWYodiA8IDApIHtcbiAgICAgICAgICBpZigrK2NvdW50ID49IDIpIHtcbiAgICAgICAgICAgIHJldHVybiBmYWxzZVxuICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBjZWxsW2pdID0gdlxuICAgICAgfVxuICAgICAgcmV0dXJuIHRydWVcbiAgICB9KVxuICB9IGVsc2Uge1xuICAgIGh1bGwgPSBodWxsLmZpbHRlcihmdW5jdGlvbihjZWxsKSB7XG4gICAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICAgIHZhciB2ID0gZGluZGV4W2NlbGxbaV1dXG4gICAgICAgIGlmKHYgPCAwKSB7XG4gICAgICAgICAgcmV0dXJuIGZhbHNlXG4gICAgICAgIH1cbiAgICAgICAgY2VsbFtpXSA9IHZcbiAgICAgIH1cbiAgICAgIHJldHVybiB0cnVlXG4gICAgfSlcbiAgfVxuXG4gIGlmKGQgJiAxKSB7XG4gICAgZm9yKHZhciBpPTA7IGk8aHVsbC5sZW5ndGg7ICsraSkge1xuICAgICAgdmFyIGggPSBodWxsW2ldXG4gICAgICB2YXIgeCA9IGhbMF1cbiAgICAgIGhbMF0gPSBoWzFdXG4gICAgICBoWzFdID0geFxuICAgIH1cbiAgfVxuXG4gIHJldHVybiBodWxsXG59IiwiXG5tb2R1bGUuZXhwb3J0cyA9IHBhcnNlXG5cbi8qKlxuICogZXhwZWN0ZWQgYXJndW1lbnQgbGVuZ3Roc1xuICogQHR5cGUge09iamVjdH1cbiAqL1xuXG52YXIgbGVuZ3RoID0ge2E6IDcsIGM6IDYsIGg6IDEsIGw6IDIsIG06IDIsIHE6IDQsIHM6IDQsIHQ6IDIsIHY6IDEsIHo6IDB9XG5cbi8qKlxuICogc2VnbWVudCBwYXR0ZXJuXG4gKiBAdHlwZSB7UmVnRXhwfVxuICovXG5cbnZhciBzZWdtZW50ID0gLyhbYXN0dnpxbWhsY10pKFteYXN0dnpxbWhsY10qKS9pZ1xuXG4vKipcbiAqIHBhcnNlIGFuIHN2ZyBwYXRoIGRhdGEgc3RyaW5nLiBHZW5lcmF0ZXMgYW4gQXJyYXlcbiAqIG9mIGNvbW1hbmRzIHdoZXJlIGVhY2ggY29tbWFuZCBpcyBhbiBBcnJheSBvZiB0aGVcbiAqIGZvcm0gYFtjb21tYW5kLCBhcmcxLCBhcmcyLCAuLi5dYFxuICpcbiAqIEBwYXJhbSB7U3RyaW5nfSBwYXRoXG4gKiBAcmV0dXJuIHtBcnJheX1cbiAqL1xuXG5mdW5jdGlvbiBwYXJzZShwYXRoKSB7XG5cdHZhciBkYXRhID0gW11cblx0cGF0aC5yZXBsYWNlKHNlZ21lbnQsIGZ1bmN0aW9uKF8sIGNvbW1hbmQsIGFyZ3Mpe1xuXHRcdHZhciB0eXBlID0gY29tbWFuZC50b0xvd2VyQ2FzZSgpXG5cdFx0YXJncyA9IHBhcnNlVmFsdWVzKGFyZ3MpXG5cblx0XHQvLyBvdmVybG9hZGVkIG1vdmVUb1xuXHRcdGlmICh0eXBlID09ICdtJyAmJiBhcmdzLmxlbmd0aCA+IDIpIHtcblx0XHRcdGRhdGEucHVzaChbY29tbWFuZF0uY29uY2F0KGFyZ3Muc3BsaWNlKDAsIDIpKSlcblx0XHRcdHR5cGUgPSAnbCdcblx0XHRcdGNvbW1hbmQgPSBjb21tYW5kID09ICdtJyA/ICdsJyA6ICdMJ1xuXHRcdH1cblxuXHRcdHdoaWxlICh0cnVlKSB7XG5cdFx0XHRpZiAoYXJncy5sZW5ndGggPT0gbGVuZ3RoW3R5cGVdKSB7XG5cdFx0XHRcdGFyZ3MudW5zaGlmdChjb21tYW5kKVxuXHRcdFx0XHRyZXR1cm4gZGF0YS5wdXNoKGFyZ3MpXG5cdFx0XHR9XG5cdFx0XHRpZiAoYXJncy5sZW5ndGggPCBsZW5ndGhbdHlwZV0pIHRocm93IG5ldyBFcnJvcignbWFsZm9ybWVkIHBhdGggZGF0YScpXG5cdFx0XHRkYXRhLnB1c2goW2NvbW1hbmRdLmNvbmNhdChhcmdzLnNwbGljZSgwLCBsZW5ndGhbdHlwZV0pKSlcblx0XHR9XG5cdH0pXG5cdHJldHVybiBkYXRhXG59XG5cbmZ1bmN0aW9uIHBhcnNlVmFsdWVzKGFyZ3Mpe1xuXHRhcmdzID0gYXJncy5tYXRjaCgvLT9bLjAtOV0rKD86ZVstK10/XFxkKyk/L2lnKVxuXHRyZXR1cm4gYXJncyA/IGFyZ3MubWFwKE51bWJlcikgOiBbXVxufVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgc2lnbiA9IHJlcXVpcmUoJy4vdXRpbGl0aWVzLmpzJykuc2lnbjtcbnZhciBjYWxjdWxhdGVEaXN0YW5jZSA9IHJlcXVpcmUoJy4vdXRpbGl0aWVzLmpzJykuZGlzdGFuY2U7XG5cbi8vIHZhciBwb2ludHMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5wb2ludHM7XG4vLyB2YXIgY2l0eVNldCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLmNpdHlTZXQ7XG4vLyB2YXIgdGV4dFBvaW50c0lkID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykudGV4dFBvaW50c0lkO1xuLy8gdmFyIHBvc3NpYmxlU3RhcnRQb2ludHNJZCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnBvc3NpYmxlU3RhcnRQb2ludHNJZDtcblxudmFyIGxpdmVNb3VzZVBvc2l0aW9uID0gcmVxdWlyZSgnLi9tb3VzZS5qcycpO1xuXG52YXIgVmVjdG9yID0gcmVxdWlyZSgnLi92ZWN0b3IuanMnKTtcblxudmFyIHJhbmRvbSA9IE1hdGgucmFuZG9tO1xudmFyIGZsb29yID0gTWF0aC5mbG9vcjtcbi8vIHZhciBSRVBVTFNJT04gPSAwLjA1O1xuLy8gdmFyIFJFUFVMU0lPTlNQRUVEID0gMC4wMDI7XG4vLyB2YXIgQU5UVkVMT0NJVFkgPSAwLjAwMTtcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihjb250YWluZXIsIGluaXRQb2ludHMsIG9wdGlvbnMpe1xuXG4gICAgY29uc29sZS5sb2coJ09wdGlvbnMgYW50IDonLCBvcHRpb25zKTtcbiAgICAvLyBEZWZpbmUgdGhvc2UgcGFyYW1ldGVycyBhcyBhdHRyaWJ1dGVzIG9mIEFudCBvYmplY3QgP1xuICAgIHZhciBSRVBVTFNJT04gPSBvcHRpb25zLnJlcFNpemU7XG4gICAgdmFyIFJFUFVMU0lPTlNQRUVEID0gb3B0aW9ucy5yZXBTcGVlZDtcbiAgICB2YXIgQU5UVkVMT0NJVFkgPSBvcHRpb25zLnZlbG9jaXR5O1xuICAgIHZhciBXRUlHSFQgPSBvcHRpb25zLndlaWdodDtcblxuICAgIHZhciBtb3VzZSA9IGxpdmVNb3VzZVBvc2l0aW9uKGNvbnRhaW5lcik7XG5cbiAgICB2YXIgcG9pbnRzID0gaW5pdFBvaW50cy5wb2ludHM7XG4gICAgdmFyIGNpdHlTZXQgPSBpbml0UG9pbnRzLmNpdHlTZXQ7XG4gICAgdmFyIHRleHRQb2ludHNJZCA9IGluaXRQb2ludHMudGV4dFBvaW50c0lkO1xuICAgIHZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSBpbml0UG9pbnRzLnBvc3NpYmxlU3RhcnRQb2ludHNJZDtcblxuXG4gICAgZnVuY3Rpb24gQW50KHBvaW50KSB7XG4gICAgICAgIHRoaXMueCA9IHBvaW50Lng7ICAgICAgICAgICAgICAgIFxuICAgICAgICB0aGlzLnkgPSBwb2ludC55O1xuICAgICAgICB0aGlzLnZlbG9jaXR5ID0gQU5UVkVMT0NJVFk7XG4gICAgICAgIHRoaXMud2VpZ2h0ID0gV0VJR0hUO1xuICAgICAgICB0aGlzLnJlcFNpemUgPSBSRVBVTFNJT047XG4gICAgICAgIHRoaXMucmVwU3BlZWQgPSBSRVBVTFNJT05TUEVFRDtcbiAgICAgICAgdGhpcy5lZGdlID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLnN0YXRlID0gXCJmb3JhZ2VcIjtcbiAgICAgICAgdGhpcy5lZGdlcyA9IFtdO1xuICAgICAgICB0aGlzLmxhc3RDaXR5ID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLm9yaWdpbiA9IHBvaW50O1xuICAgICAgICB0aGlzLmRlc3RpbmF0aW9uID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLm9yaWVudGF0aW9uID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLmRpcmVjdGlvbiA9IG5ldyBWZWN0b3IoMCwwKTtcbiAgICAgICAgdGhpcy5wcm9nID0gMDtcbiAgICB9XG4gICAgLy8gZm9yYWdlOiB0aGUgYW50IHdhbmRlcnMgYXJvdW5kIHdpdGhvdXQgYW55IHBoZXJvbW9uIGRlcG9zaXRpb25cbiAgICAvLyBvbmNlIGl0IGZpbmRzIGEgY2l0eSwgaXQgc3RhcnRzIHJlbWVtYmVyaW5nIHRoZSBub2RlcyBpdCBnb2VzIHRocm91Z2hcbiAgICAvLyB3aGVuIGl0IGZpbmRzIGFub3RoZXIgY2l0eSwgaXQgY29tcHV0ZXMgdGhlIHBhdGggbGVuZ3RoIGFuZCBhZGRzIHBoZXJvbW9ucyBvbmUgZWFjaCBlZGdlc1xuICAgIC8vIHByb3BvcnRpb25uYWx5IHRvIHRoZSBzaG9ydGVzdG5lc3Mgb2YgdGhlIHBhdGhcbiAgICAvLyBpdCByZXNldHMgdGhlIGxpc3Qgb2Ygbm9kZXMgYW5kIGNvbnRpbnVlc1xuICAgIC8vIHdoaWxlIGZvcmFnaW5nIHRoZSBhbnQgY2hvc2VzIHRoZSBwYXRoIHdpdGggYSBwaGVyb21vbiBwcmVmZXJlbmNlXG5cblxuICAgIC8vIHN0YXRpYyBtZXRob2RzXG4gICAgQW50LmdlbmVyYXRlUmFuZFN0YXJ0UG9pbnQgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgdmFyIHJhbmRJZCA9IE1hdGguZmxvb3IocG9zc2libGVTdGFydFBvaW50c0lkLmxlbmd0aCAqIHJhbmRvbSgpKTtcbiAgICAgICAgdmFyIHJhbmRTdGFydFBvaW50ID0gcG9pbnRzW3Bvc3NpYmxlU3RhcnRQb2ludHNJZFtyYW5kSWRdXTtcbiAgICAgICAgcmV0dXJuIHJhbmRTdGFydFBvaW50O1xuICAgIH1cblxuXG4gICAgLy8gbWV0aG9kc1xuICAgIEFudC5wcm90b3R5cGUgPSB7XG5cbiAgICAgICAgdHJhbnNpdDogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHN3aXRjaCAodGhpcy5zdGF0ZSkge1xuICAgICAgICAgICAgY2FzZSBcImZvcmFnZVwiOlxuICAgICAgICAgICAgICAgIHZhciByZXMgPSB0aGlzLm1vdmUoKTtcbiAgICAgICAgICAgICAgICBpZiAocmVzLmNpdHlSZWFjaGVkKSB7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMuc3RhdGUgPSBcInBoZXJvbW9uaW5nXCI7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMubGFzdENpdHkgPSB0aGlzLm9yaWdpbi5pZDtcbiAgICAgICAgICAgICAgICB9O1xuICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgY2FzZSBcInBoZXJvbW9uaW5nXCI6XG4gICAgICAgICAgICAgICAgdmFyIHJlcyA9IHRoaXMubW92ZSgpO1xuICAgICAgICAgICAgICAgIGlmIChyZXMuZWRnZUNoYW5nZWQpIHtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5lZGdlcy5wdXNoKHRoaXMuZWRnZSk7XG4gICAgICAgICAgICAgICAgICAgIC8vIGZvdW5kIGEgY2l0eVxuICAgICAgICAgICAgICAgICAgICBpZiAocmVzLmNpdHlSZWFjaGVkICYmICh0aGlzLm9yaWdpbi5pZCAhPSB0aGlzLmxhc3RDaXR5KSApe1xuICAgICAgICAgICAgICAgICAgICAgICAgLy8gY29tcHV0ZSB0aGUgbGVuZ3RoIG9mIHRoZSBwYXRoXG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgcGF0aExlbmd0aCA9IHRoaXMuZWRnZXMubWFwKGZ1bmN0aW9uKGUpe3JldHVybiBlLmRpc3RhbmNlfSkucmVkdWNlKGZ1bmN0aW9uKGEsYil7cmV0dXJuIGEgKyBifSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgZGVsdGFQaGVyb21vbmUgPSAxL3BhdGhMZW5ndGg7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgYW50V2VpZ2h0ID0gdGhpcy53ZWlnaHQ7XG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzLmZvckVhY2goZnVuY3Rpb24oZSl7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGEgPSBlLnB0MSwgYiA9IGUucHQyLCB3ZWlnaHQgPSAxOyAgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgLy8gaW5jcmVhc2VkIGRyb3BwZWQgcGhlcm9tb25zIGZvciB0ZXh0RWRnZXNcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoKGNpdHlTZXQuaW5kZXhPZihhLmlkKSAhPSAtMSkgJiYgY2l0eVNldC5pbmRleE9mKGIuaWQpICE9IC0xICYmIChNYXRoLmFicyhhLmlkIC0gYi5pZCkgPT0gMSkpXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB3ZWlnaHQgKj0gYW50V2VpZ2h0O1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBlLnBoZXJvbW9uICs9IChkZWx0YVBoZXJvbW9uZSAqIHdlaWdodCk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9KTtcblxuICAgICAgICAgICAgICAgICAgICAgICAgdGhpcy5lZGdlcyA9IFt0aGlzLmVkZ2VdO1xuICAgICAgICAgICAgICAgICAgICAgICAgdGhpcy5sYXN0Q2l0eSA9IHRoaXMub3JpZ2luLmlkO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgIH1cblxuICAgICAgICB9LFxuXG4gICAgICAgIHNldERpcmVjdGlvbjogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHZhciBwb3NzaWJsZUVkZ2VzID0gW107XG5cbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdGhpcy5vcmlnaW4ubmV4dHMubGVuZ3RoOyBpKyspXG4gICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgcG9zc2libGVFZGdlc1tpXSA9IHRoaXMub3JpZ2luLm5leHRzW2ldO1xuICAgICAgICAgICAgfSBcblxuICAgICAgICAgICAgcG9zc2libGVFZGdlcy5zcGxpY2UocG9zc2libGVFZGdlcy5pbmRleE9mKHRoaXMuZWRnZSksMSk7XG5cbiAgICAgICAgICAgIC8vIGZsaXAgYSBjb2luIGFuZCBlaXRoZXIgdGFrZSB0aGUgc21lbGxpZXN0IHBhdGggb3IgYSByYW5kb20gb25lXG4gICAgICAgICAgICBpZiAocmFuZG9tKCkgPiAwLjUpe1xuICAgICAgICAgICAgICAgIHZhciBzbWVsbHMgPSBwb3NzaWJsZUVkZ2VzLm1hcChmdW5jdGlvbihlKXtyZXR1cm4gZS5waGVyb21vbn0pO1xuICAgICAgICAgICAgICAgIHZhciBpbmRleCA9IHNtZWxscy5pbmRleE9mKE1hdGgubWF4LmFwcGx5KE1hdGgsIHNtZWxscykpO1xuICAgICAgICAgICAgICAgIHRoaXMuZWRnZSA9IHBvc3NpYmxlRWRnZXNbaW5kZXhdO1xuICAgICAgICAgICAgfSBcbiAgICAgICAgICAgIGVsc2VcbiAgICAgICAgICAgICAgICB0aGlzLmVkZ2UgPSBwb3NzaWJsZUVkZ2VzW2Zsb29yKHJhbmRvbSgpKnBvc3NpYmxlRWRnZXMubGVuZ3RoKV07XG5cbiAgICAgICAgICAgIC8vIHNldCB0aGUgZGVzdGluYXRpb24gcG9pbnQsIGJlaW5nIGVkZ2UucHQxIG9yIGVkZ2UucHQyXG4gICAgICAgICAgICB0aGlzLmRlc3RpbmF0aW9uID0gKHRoaXMub3JpZ2luID09IHRoaXMuZWRnZS5wdDEpID8gdGhpcy5lZGdlLnB0MiA6IHRoaXMuZWRnZS5wdDE7XG5cbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLnggPSB0aGlzLmRlc3RpbmF0aW9uLnggLSB0aGlzLm9yaWdpbi54OyBcbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLnkgPSB0aGlzLmRlc3RpbmF0aW9uLnkgLSB0aGlzLm9yaWdpbi55O1xuXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi5ub3JtYWxpemUoKTtcbiAgICAgICAgfSxcblxuICAgICAgICBtb3ZlOiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgdmFyIGVkZ2VDaGFuZ2VkO1xuICAgICAgICAgICAgdmFyIGNpdHlSZWFjaGVkID0gZmFsc2U7XG5cbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLnggPSB0aGlzLmRlc3RpbmF0aW9uLnggLSB0aGlzLng7IFxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueSA9IHRoaXMuZGVzdGluYXRpb24ueSAtIHRoaXMueTtcbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLm5vcm1hbGl6ZSgpO1xuXG4gICAgICAgICAgICAvLyBvbiBlZGdlXG4gICAgICAgICAgICBpZiAoKGNhbGN1bGF0ZURpc3RhbmNlKHRoaXMsIHRoaXMuZGVzdGluYXRpb24pID4gdGhpcy5yZXBTcGVlZCkpe1xuXG4gICAgICAgICAgICAgICAgLy8gYSBkZWx0YSBtb3ZlbWVudCB3aWxsIGJlIGFwcGxpZWQgaWYgY29sbGlzaW9uIHdpdGggb2JzdGFjbGUgZGV0ZWN0ZWRcbiAgICAgICAgICAgICAgICB2YXIgZGVsdGEgPSB0aGlzLmF2b2lkT2JzdGFjbGUoKTtcblxuICAgICAgICAgICAgICAgIHRoaXMueCArPSB0aGlzLnZlbG9jaXR5ICogdGhpcy5kaXJlY3Rpb24ueCArIGRlbHRhLnggKiB0aGlzLnJlcFNwZWVkO1xuICAgICAgICAgICAgICAgIHRoaXMueSArPSB0aGlzLnZlbG9jaXR5ICogdGhpcy5kaXJlY3Rpb24ueSArIGRlbHRhLnkgKiB0aGlzLnJlcFNwZWVkO1xuXG4gICAgICAgICAgICAgICAgdGhpcy5wcm9nID0gdGhpcy5jYWxjdWxhdGVQcm9ncmVzc2lvbigpO1xuICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgIGVkZ2VDaGFuZ2VkID0gZmFsc2U7XG5cbiAgICAgICAgICAgIC8vIG9uIHZlcnRleFxuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICB0aGlzLnN0ZXAgPSAwO1xuICAgICAgICAgICAgICAgIHRoaXMucHJvZyA9IDA7XG4gICAgICAgICAgICAgICAgdGhpcy5vcmlnaW4gPSB0aGlzLmRlc3RpbmF0aW9uO1xuICAgICAgICAgICAgICAgIHRoaXMueCA9IHRoaXMub3JpZ2luLng7XG4gICAgICAgICAgICAgICAgdGhpcy55ID0gdGhpcy5vcmlnaW4ueTtcblxuICAgICAgICAgICAgICAgIHRoaXMuc2V0RGlyZWN0aW9uKCk7XG5cbiAgICAgICAgICAgICAgICBjaXR5UmVhY2hlZCA9IChjaXR5U2V0LmluZGV4T2YodGhpcy5vcmlnaW4uaWQpICE9IC0xKTtcbiAgICAgICAgICAgICAgICBlZGdlQ2hhbmdlZCA9IHRydWU7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXR1cm4ge2NpdHlSZWFjaGVkOiBjaXR5UmVhY2hlZCwgZWRnZUNoYW5nZWQ6IGVkZ2VDaGFuZ2VkfTtcbiAgICAgICAgfSxcblxuICAgICAgICBhdm9pZE9ic3RhY2xlOiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgdmFyIGRpc3RhbmNlID0gY2FsY3VsYXRlRGlzdGFuY2UodGhpcywgbW91c2UpO1xuICAgICAgICBcbiAgICAgICAgICAgIGlmIChkaXN0YW5jZSA8PSB0aGlzLnJlcFNpemUpIHtcblxuICAgICAgICAgICAgICAgIHJldHVybiB7XG4gICAgICAgICAgICAgICAgICAgIC8vIGRlbHRhIG1vdmVtZW50IGlzIGNvbXBvc2VkIG9mIGEgcmVwdWxzaW9uIGRlbHRhIGFuZCBhIGNpcmN1bGFyIGRlbHRhIFxuICAgICAgICAgICAgICAgICAgICB4OiAodGhpcy54IC0gbW91c2UueCkvZGlzdGFuY2UgKyAodGhpcy55IC0gbW91c2UueSkvZGlzdGFuY2UgKiAxLFxuICAgICAgICAgICAgICAgICAgICB5OiAodGhpcy55IC0gbW91c2UueSkvZGlzdGFuY2UgLSAodGhpcy54IC0gbW91c2UueCkvZGlzdGFuY2UgKiAxXG4gICAgICAgICAgICAgICAgfTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGVsc2VcbiAgICAgICAgICAgICAgICByZXR1cm4ge3g6MCwgeTowfTtcbiAgICAgICAgfSxcblxuICAgICAgICBjYWxjdWxhdGVQcm9ncmVzc2lvbjogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHZhciB2ID0gbmV3IFZlY3Rvcih0aGlzLnggLSB0aGlzLm9yaWdpbi54LCB0aGlzLnkgLSB0aGlzLm9yaWdpbi55KTtcbiAgICAgICAgICAgIHZhciBub3JtID0gdi5ub3JtKCk7XG5cbiAgICAgICAgICAgIHZhciB0aGV0YSA9ICh2LnggKiB0aGlzLmVkZ2UuZGlyZWN0aW9uLnggKyB2LnkgKiB0aGlzLmVkZ2UuZGlyZWN0aW9uLnkpIC8gbm9ybTtcbiAgICAgICAgICAgIHZhciBwcm9nID0gbm9ybSAqIE1hdGguYWJzKHRoZXRhKTtcbiAgICAgICAgICAgIC8vIHJldHVybnMgbGVuZ3RoIG9mIHByb2plY3Rpb24gb24gZWRnZVxuICAgICAgICAgICAgcmV0dXJuIHByb2c7XG4gICAgICAgIH1cblxuICAgIH07XG4gICAgcmV0dXJuIEFudDtcbn1cblxuIiwiJ3VzZSBzdHJpY3QnXG5cbnZhciBkdCA9IHJlcXVpcmUoXCJkZWxhdW5heS10cmlhbmd1bGF0ZVwiKTtcblxudmFyIHJhbmdlID0gcmVxdWlyZSgnLi91dGlsaXRpZXMuanMnKS5yYW5nZTtcblxuLy8gdmFyIHBvaW50cyA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnBvaW50cztcbnZhciB0ZXh0TWVzaCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnRleHRNZXNoO1xudmFyIGNpdHlTZXQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5jaXR5U2V0O1xudmFyIG5iUmFuZG9tUG9pbnRzID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykubmJSYW5kb21Qb2ludHM7XG52YXIgZm9yY2VkRWRnZXMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5mb3JjZWRFZGdlcztcblxudmFyIEVkZ2UgPSByZXF1aXJlKCcuL2VkZ2UuanMnKTtcblxuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKHBvaW50cyl7XG4gICAgLy8gdHJpYW5ndWxhdGVcbiAgICB2YXIgY2VsbHMgPSBkdChwb2ludHMubWFwKGZ1bmN0aW9uKHApe1xuICAgICAgICByZXR1cm4gW3AueCwgcC55XTtcbiAgICB9KSk7XG5cbiAgICB2YXIgZWRnZXMgPSBbXTtcbiAgICB2YXIgcGVybXV0YXRpb25zID0gW1swLDFdLCBbMCwyXSwgWzEsMl1dO1xuXG4gICAgLy8gZm9yY2UgdGhlIGVkZ2VzIG9mIHRoZSB0ZXh0IHRvIGJlIGVkZ2VzIG9mIHRoZSBncmFwaFxuICAgIGlmICh0ZXh0TWVzaCkge1xuICAgICAgICByYW5nZSgwLCBwb2ludHMubGVuZ3RoIC0gbmJSYW5kb21Qb2ludHMpLmZvckVhY2goZnVuY3Rpb24oaWQpe1xuICAgICAgICAgICAgdmFyIGRpcmVjdExpbmsgPSBmb3JjZWRFZGdlc1tpZF07XG4gICAgICAgICAgICB2YXIgdGV4dEVkZ2UgPSBFZGdlLmNyZWF0ZShwb2ludHNbaWRdLCBwb2ludHNbZGlyZWN0TGlua10pO1xuICAgICAgICAgICAgZWRnZXMucHVzaCh0ZXh0RWRnZSk7XG4gICAgICAgICAgICBwb2ludHNbaWRdLm5leHRzLnB1c2godGV4dEVkZ2UpO1xuICAgICAgICB9KVxuICAgIH1cblxuXG4gICAgY2VsbHMuZm9yRWFjaChmdW5jdGlvbihjZWxsKXtcbiAgICAgICBcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCAzOyArK2kpeyAgLy8gZm9yIGVhY2ggcG9pbnQuaWQgbGlzdGVkIGluIGN1cnJlbnQgY2VsbFxuICAgICAgICAgICAgdmFyIHB0ID0gcG9pbnRzW2NlbGxbaV1dO1xuXG4gICAgICAgICAgICBmb3IgKHZhciBqID0gMTsgaiA8IDM7ICsrail7IFxuXG4gICAgICAgICAgICAgICAgdmFyIHB0aiA9IHBvaW50c1tjZWxsWyggaSArIGogKSAlIDNdXTsgLy8gcGljayBvbmUgb2YgdGhlIG90aGVyIDIgcG9pbnRzIG9mIHRoZSBjZWxsXG4gICAgICAgICAgICAgICAgdmFyIG5ld0VkZ2UgPSB1bmRlZmluZWQ7XG5cbiAgICAgICAgICAgICAgICAvLyBpZiBwdCBhbHJlYWR5IGhhcyBuZXh0RWRnZXMgLi4uXG4gICAgICAgICAgICAgICAgaWYgKHB0Lm5leHRzLmxlbmd0aCAhPSAwKSB7XG4gICAgICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgICAgICAvLyAuLi4gZ2V0IHRoZSBwb2ludHMgY29ycmVzcG9uZGluZyAuLi5cbiAgICAgICAgICAgICAgICAgICAgdmFyIHRlbXBQb2ludHMgPSBwdC5uZXh0cy5tYXAoZnVuY3Rpb24oZSl7XG4gICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gW2UucHQxLCBlLnB0Ml07XG4gICAgICAgICAgICAgICAgICAgIH0pLnJlZHVjZShmdW5jdGlvbihhLCBiKXtcbiAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gYS5jb25jYXQoYik7XG4gICAgICAgICAgICAgICAgICAgIH0pO1xuXG4gICAgICAgICAgICAgICAgICAgIC8vIC4uLiBhbmQgY2hlY2sgaWYgcHRqIGFscmVhZHkgaXMgcGFydCBvZiB0aGUgZXhpc3RpbmcgbmV4dEVkZ2VzLiBJZiBub3QsIGFkZCB0aGUgZWRnZS5cbiAgICAgICAgICAgICAgICAgICAgaWYgKHRlbXBQb2ludHMuaW5kZXhPZihwdGopID09IC0xKXtcbiAgICAgICAgICAgICAgICAgICAgICAgIG5ld0VkZ2UgPSBFZGdlLmNyZWF0ZShwdCwgcHRqKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGVkZ2VzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgICAgICAgICBwdC5uZXh0cy5wdXNoKG5ld0VkZ2UpO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICBuZXdFZGdlID0gRWRnZS5jcmVhdGUocHQsIHB0aik7XG4gICAgICAgICAgICAgICAgICAgIGVkZ2VzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgICAgIHB0Lm5leHRzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAgICAgLy8gYWRkIGFsc28gdGhlIGVkZ2UgdG8gdGhlIGVkZ2UncyBvdGhlciBwb2ludCdzIG5leHRFZGdlc1xuICAgICAgICAgICAgICAgIGlmIChuZXdFZGdlICE9IHVuZGVmaW5lZCl7XG4gICAgICAgICAgICAgICAgICAgIHB0ai5uZXh0cy5wdXNoKG5ld0VkZ2UpO1xuICAgICAgICAgICAgICAgIH0gICAgICAgICBcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgLy8gYWRkIHRoZSB0ZXh0RWRnZXMgdG8gbmV4dEVkZ2VzIG1hcFxuICAgICAgICAgICAgaWYgKHRleHRNZXNoICYmIChjaXR5U2V0LmluZGV4T2YocHQpICE9IC0xKSkge1xuICAgICAgICAgICAgICAgIHZhciB0ZXh0RWRnZSA9IEVkZ2UuY3JlYXRlKHB0LCBwb2ludHNbcHQuaWQgKyAxXSk7XG4gICAgICAgICAgICAgICAgZWRnZXMucHVzaCh0ZXh0RWRnZSk7XG4gICAgICAgICAgICAgICAgcHQubmV4dHMucHVzaCh0ZXh0RWRnZSk7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgfVxuICAgIH0pO1xuXG4gICAgcmV0dXJuIGVkZ2VzO1xufTsiLCIndXNlIHN0cmljdCc7XG5cbnZhciBzcXJ0ID0gTWF0aC5zcXJ0O1xudmFyIHBvdyA9IE1hdGgucG93O1xudmFyIGFicyA9IE1hdGguYWJzO1xudmFyIGF0YW4gPSBNYXRoLmF0YW47XG5cbnZhciBWZWN0b3IgPSByZXF1aXJlKCcuL3ZlY3Rvci5qcycpO1xuXG5cbmZ1bmN0aW9uIEVkZ2UocHRBLCBwdEIpIHtcbiAgICB2YXIgZGlzdGFuY2UgPSBzcXJ0KCBwb3cocHRBLnggLSBwdEIueCwgMikgKyBwb3cocHRBLnkgLSBwdEIueSwgMikgKTtcblxuICAgIC8vIGZpbmQgbGluZSBlcXVhdGlvbiBheCArIGJ5ICsgYyA9IDBcbiAgICB2YXIgYSA9IDE7XG4gICAgdmFyIGIgPSAtIChwdEIueCAtIHB0QS54KSAvIChwdEIueSAtIHB0QS55KTtcblxuICAgIC8vIG9yaWVudGF0ZSB2ZWN0b3IgKGEsYilcbiAgICBpZiAoYiA8IDApe1xuICAgICAgICBiID0gLWI7XG4gICAgICAgIGEgPSAtYTtcbiAgICB9XG5cbiAgICAvLyBub3JtYWxpemUgdmVjdG9yIChhLGIpXG4gICAgdmFyIG4gPSBuZXcgVmVjdG9yKGEsIGIpO1xuICAgIG4ubm9ybWFsaXplKCk7XG5cbiAgICB2YXIgYyA9IC0gKGEgKiBwdEEueCArIGIgKiBwdEEueSk7XG5cbiAgICAvLyAvLyBjYWxjdWxhdGUgdmVjdG9yIGRpcmVjdG9yXG4gICAgdmFyIHYgPSBuZXcgVmVjdG9yKHB0Qi54IC0gcHRBLngsIHB0Qi55IC0gcHRBLnkpO1xuICAgIFxuICAgIHYubm9ybWFsaXplKCk7XG5cblxuICAgIHRoaXMuaWQgPSB1bmRlZmluZWQ7XG4gICAgdGhpcy5wdDEgPSBwdEE7XG4gICAgdGhpcy5wdDIgPSBwdEI7XG4gICAgdGhpcy5kaXJlY3Rpb24gPSB2O1xuICAgIHRoaXMub3J0aERpcmVjdGlvbiA9IG47IFxuICAgIHRoaXMuZGlzdGFuY2UgPSBkaXN0YW5jZTtcbiAgICB0aGlzLnBoZXJvbW9uID0gMS9kaXN0YW5jZTtcbiAgICB0aGlzLmxpbmUgPSB7XG4gICAgICAgIGE6IGEsXG4gICAgICAgIGI6IGIsXG4gICAgICAgIGM6IGMsXG4gICAgfVxufVxuXG5cbi8vIHN0YXRpYyBtZXRob2RzXG5FZGdlLmNyZWF0ZSA9IGZ1bmN0aW9uKHB0QSwgcHRCKSB7XG4gICAgdmFyIGVkZ2UgPSBuZXcgRWRnZShwdEEsIHB0Qik7XG4gICAgcmV0dXJuIGVkZ2U7XG59XG5cblxuLy8gbWV0aG9kc1xuRWRnZS5wcm90b3R5cGUgPSB7XG5cbiAgICBnZXRPdGhlclBvaW50OiBmdW5jdGlvbihwb2ludCkge1xuICAgICAgICBpZiAocG9pbnQgPT0gdGhpcy5wdDEpXG4gICAgICAgICAgICByZXR1cm4gdGhpcy5wdDI7XG4gICAgICAgIGVsc2UgaWYgKHBvaW50ID09IHRoaXMucHQyKVxuICAgICAgICAgICAgcmV0dXJuIHRoaXMucHQxO1xuICAgICAgICBlbHNlXG4gICAgICAgICAgICBjb25zb2xlLmxvZyhcIkVycm9yXCIpO1xuICAgIH0sXG5cbiAgICBjYWxjdWxhdGVEaXN0YW5jZTogZnVuY3Rpb24oeCwgeSkge1xuICAgICAgICB2YXIgYSA9IHRoaXMubGluZS5hLFxuICAgICAgICAgICAgYiA9IHRoaXMubGluZS5iLFxuICAgICAgICAgICAgYyA9IHRoaXMubGluZS5jO1xuICAgICAgICByZXR1cm4gYWJzKGEgKiB4ICsgYiAqIHkgKyBjKSAvIE1hdGguc3FydChNYXRoLnBvdyhhLDIpICsgTWF0aC5wb3coYiwyKSk7XG4gICAgfSxcblxufVxubW9kdWxlLmV4cG9ydHMgPSBFZGdlOyIsIid1c2Ugc3RyaWN0J1xuXG52YXIgYW50RnVuY3Rpb24gPSByZXF1aXJlKCcuL2FudC5qcycpO1xuXG4vLyB2YXIgTkJBTlRTID0gNDAwMDtcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbiAoY29udGFpbmVyLCBwb2ludHNJbmZvcywgb3B0aW9ucykge1xuXG5cdHZhciBBbnQgPSBhbnRGdW5jdGlvbihjb250YWluZXIsIHBvaW50c0luZm9zLCBvcHRpb25zKTtcblxuXHR2YXIgcG9wdWxhdGlvbiA9IG5ldyBBcnJheShvcHRpb25zLm5iQW50cyk7XG5cdHZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSBwb2ludHNJbmZvcy5wb3NzaWJsZVN0YXJ0UG9pbnRzSWQ7XG5cblx0Zm9yICh2YXIgaSA9IDA7IGkgPCBvcHRpb25zLm5iQW50czsgaSsrKSB7XG5cdFx0dmFyIG5ld0FudCA9IG5ldyBBbnQoQW50LmdlbmVyYXRlUmFuZFN0YXJ0UG9pbnQoKSk7XG5cdFx0bmV3QW50LnNldERpcmVjdGlvbigpO1xuXHRcdHBvcHVsYXRpb25baV0gPSBuZXdBbnQ7XG5cdH1cblxuXHRyZXR1cm4gcG9wdWxhdGlvbjtcblxufSIsIid1c2Ugc3RyaWN0J1xuXG52YXIgcGFyc2UgPSByZXF1aXJlKCdwYXJzZS1zdmctcGF0aCcpO1xuXG52YXIgcmFuZ2UgPSByZXF1aXJlKCcuL3V0aWxpdGllcy5qcycpLnJhbmdlO1xuXG52YXIgUG9pbnQgPSByZXF1aXJlKCcuL3BvaW50LmpzJyk7XG5cbnZhciByYW5kb20gPSBNYXRoLnJhbmRvbTtcblxuLy8gdmFyIG5iUmFuZG9tUG9pbnRzID0gNTAwO1xuLy8gdmFyIG5iU3RhcnRQb2ludHMgPSAyMDtcblxudmFyIG5iQ2l0eSA9IDI7XG5cbnZhciB0ZXh0TWVzaCA9IHRydWU7XG5cbi8vIEZyYW1lIGRlZmluaXRpb25cbnZhciB4SW5pdCA9IDAsIHlJbml0ID0gMDtcbnZhciB3ID0gMSxcbiAgICBoID0gMTtcblxudmFyIEFjaGFyID0gXCJjIDQuNjAxMSwtMTEuNzEwNDcgOS4yMDgzNSwtMjMuNDIwMDYgMTMuODE5OSwtMzUuMTI4OTggNC42MTE1NiwtMTEuNzA4OTIgOS4yMjc0MSwtMjMuNDE3MTggMTMuODQ1NzMsLTM1LjEyNSA0LjYxODMxLC0xMS43MDc4MiA5LjIzOTA4LC0yMy40MTUxOSAxMy44NjA0NiwtMzUuMTIyMzMgNC42MjEzOCwtMTEuNzA3MTQgOS4yNDMzNiwtMjMuNDE0MDYgMTMuODY0MTEsLTM1LjEyMDk3IDQuNjIwNzQsLTExLjcwNjkxIDkuMjQwMjUsLTIzLjQxMzgxIDEzLjg1NjY3LC0zNS4xMjA5MiA0LjYxNjQxLC0xMS43MDcxMiA5LjIyOTc0LC0yMy40MTQ0NCAxMy44MzgxNCwtMzUuMTIyMTggNC42MDgzOSwtMTEuNzA3NzUgOS4yMTE4NiwtMjMuNDE1OTIgMTMuODA4NTQsLTM1LjEyNDc0IDQuNTk2NjgsLTExLjcwODgxIDkuMTg2NTgsLTIzLjQxODI3NyAxMy43Njc4NSwtMzUuMTI4NjAzIDEyLjk5OTIzLC0zLjM1Nzg1NSAyNC4zMDA2OSwtNi44MjExNDQgMzMuODY4OTMsLTQuNDY0MTEgOS41NjgyNSwyLjM1NzAzNSAxNy40MDMyOCwxMC41MzQzOTMgMjMuNDY5NjcsMzAuNDU3ODMzIDQuMTc1OTgsMTAuNjIxMTQgOC4zNjAzMSwyMS4yMzg4MiAxMi41NDk0MiwzMS44NTQ1MyA0LjE4OTEsMTAuNjE1NzEgOC4zODI5OCwyMS4yMjk0MyAxMi41NzgwNywzMS44NDI2NiA0LjE5NTA5LDEwLjYxMzIzIDguMzkxNCwyMS4yMjU5NSAxMi41ODUzNSwzMS44Mzk2NSA0LjE5Mzk1LDEwLjYxMzY5IDguMzg1NTUsMjEuMjI4MzYgMTIuNTcxMjMsMzEuODQ1NDkgNC4xODU2OCwxMC42MTcxMiA4LjM2NTQ0LDIxLjIzNjcgMTIuNTM1NzIsMzEuODYwMjEgNC4xNzAyOCwxMC42MjM1IDguMzMxMDgsMjEuMjUwOTMgMTIuNDc4ODMsMzEuODgzNzcgNC4xNDc3NSwxMC42MzI4MyA4LjI4MjQ1LDIxLjI3MTA3IDEyLjQwMDUzLDMxLjkxNjIgNC4xMTgwOCwxMC42NDUxMiA4LjIxOTU1LDIxLjI5NzEyIDEyLjMwMDg1LDMxLjk1NzQ5IC01LjkwODk2LDYuOTU1NjEgLTI0LjYxNjE3LDEuMTEyOTggLTM1LjU5MzcyLDMgLTE2Ljc1MjQ4LDIuODQ4NTkgLTI2Ljk2NDIxLC0wLjQxNDE2IC0yOC40MDYyOCwtMTkgLTMuMTc3MjYsLTguNDIxNTcgLTYuMzU2MDYsLTE2Ljg3OTI0IC05LjQ3OTk3LC0yNS4zNDg1NCAtMy4xMjM5MSwtOC40NjkyOSAtNi4xOTI5MywtMTYuOTUwMiAtOS4xNTA2MiwtMjUuNDE4MjYgLTE2LjQ2NzIzLC0wLjgwMDUzIC0zNS4xNzc2OCwtMi43NDI4MSAtNTMuMzM3MjcsLTMuMTE0MzkgLTE4LjE1OTYsLTAuMzcxNTggLTM1Ljc2ODM0LDAuODI3NTMgLTUwLjAzMjE0LDYuMzA5NzYgLTQuMTU2NzksMTAuOTMxMjQgLTguMTM4MDYsMjEuOTI2OSAtMTIuMTU3ODUsMzIuOTA1NDYgLTQuMDE5OCwxMC45Nzg1NSAtOC4wNzgxMywyMS45NDAwMSAtMTIuMzg5MDMsMzIuODAyODIgLTkuNTI3NTgsMC44Mjc0NyAtMTkuMTA0NTcsMC45NjQyMyAtMjguNjkyODIsMC45MzM2MyAtOS41ODgyNCwtMC4wMzA2IC0xOS4xODc3MywtMC4yMjg1NSAtMjguNzYwMywtMC4wNzA1IGwgMCwtMi4xOTc5MSB6IG0gMTc0LC0xMDkgYyAtMy4xNzQ2MiwtOS41MTEzNiAtNi40NDgzNSwtMTguOTkxNzQgLTkuNzU0MywtMjguNDYyMjQgLTYuODE1NTksLTE5LjUxNDEyIC02LjI0MjAyLC0yNS4zNzk0NiAtMTIuMzE3NCwtNDMuNjA3ODggLTMuMTY3MjIsLTkuNTE1NDMgLTEzLjYwMjU3LC0zMi4zMjg2NiAtMTYuNDk5NzMsLTQxLjkyOTg4IC02LjA4OTg5LDEyLjg4NTc2IC0xMS4xMzIsMjYuNTkyNzIgLTE1LjkzNDE1LDQwLjQ3NDc2IC00LjgwMjE1LDEzLjg4MjAzIC05LjM2NDM1LDI3LjkzOTE1IC0xNC40OTQ0Miw0MS41MjUyNCAtMi4zOTQyNSwxMi4wNTgwMSAtMjIuODgxNSw0MC4xMTExNiAyLjE4NzQ1LDM0IDExLjA1MjU2LC0wLjQ0NDg2IDIyLjUyOTM5LDAuMTEwNjEgMzMuODU2MjQsMC4yNDk1NSAxMS4zMjY4NCwwLjEzODk1IDIyLjUwMzY5LC0wLjEzODYzIDMyLjk1NjMxLC0yLjI0OTU1IHogXCI7XG52YXIgTmNoYXIgPSBcImMgMCwwIDAsLTIzLjgzMzM0IDAsLTM1Ljc1IDAsLTExLjkxNjY3IDAsLTIzLjgzMzM0IDAsLTM1Ljc1IDAsLTExLjkxNjY2IDAsLTIzLjgzMzMzIDAsLTM1Ljc1IDAsLTExLjkxNjY2NyAwLC0yMy44MzMzMzUgMCwtMzUuNzUwMDAzIDE0LjY0NzI5LDEuMDUzMTMgMzEuMDMxOTMsLTIuMDg4MDggNDQuNjA2NzksMS41NTQxNSA2LjE4MDc2LDcuNjEwOTk2IDEyLjM1NDM4LDE1LjIzMTM0MiAxOC41MTk3MywyMi44NjEwMTMgNi4xNjUzNSw3LjYyOTY3IDEyLjMyMjQ0LDE1LjI2ODY2IDE4LjQ3MDE0LDIyLjkxNjk2IDYuMTQ3Nyw3LjY0ODI5IDEyLjI4NjAyLDE1LjMwNTg4IDE4LjQxMzgzLDIyLjk3Mjc0IDYuMTI3ODEsNy42NjY4NiAxMi4yNDUxMiwxNS4zNDMgMTguMzUwODEsMjMuMDI4MzggNi4xMDU2OCw3LjY4NTM4IDEyLjE5OTc1LDE1LjM4MDAxIDE4LjI4MTA3LDIzLjA4Mzg2IDYuMDgxMzIsNy43MDM4NCAxMi4xNDk5LDE1LjQxNjkxIDE4LjIwNDYyLDIzLjEzOTE4IDYuMDU0NzIsNy43MjIyNiAxMi4wOTU1NywxNS40NTM3MiAxOC4xMjE0NSwyMy4xOTQzNSA2LjAyNTg3LDcuNzQwNjMgMTIuMDM2NzYsMTUuNDkwNDMgMTguMDMxNTYsMjMuMjQ5MzcgMC40NTU4LC0xNS40OTE0IDAuNzI2NTUsLTMwLjk4NzI1IDAuODgxMTksLTQ2LjQ4NTgyIDAuMTU0NjQsLTE1LjQ5ODU4IDAuMTkzMTksLTMwLjk5OTg4IDAuMTg0NTgsLTQ2LjUwMjE4IC0wLjAwOSwtMTUuNTAyMyAtMC4wNjQzLC0zMS4wMDU2IC0wLjA5ODMsLTQ2LjUwODE3IC0wLjAzNCwtMTUuNTAyNTggLTAuMDQ2MSwtMzEuMDA0NDMyIDAuMDMyNSwtNDYuNTAzODMzIDkuMzMzMzMsMCAxOC42NjY2NywwIDI4LDAgOS4zMzMzMywwIDE4LjY2NjY2LDAgMjgsMCAwLDExLjkxNjY2NyAwLDIzLjgzMzMzMiAwLDM1Ljc1MDAwMyAwLDExLjkxNjY2IDAsMjMuODMzMzMgMCwzNS43NSAwLDExLjkxNjY2IDAsMjMuODMzMzMgMCwzNS43NSAwLDExLjkxNjY3IDAsMjMuODMzMzMgMCwzNS43NSAwLDExLjkxNjY2IDAsMjMuODMzMzMgMCwzNS43NSAwLDExLjkxNjY3IDAsMjMuODMzMzQgMCwzNS43NSAwLDExLjkxNjY2IDAsMjMuODMzMzMgMCwzNS43NSAwLDExLjkxNjY3IDAsMjMuODMzMzQgMCwzNS43NSAtMTQuNjQ0MSwtMS4wNTM0OCAtMzEuMDI2LDIuMDg4NDcgLTQ0LjU5NzY0LC0xLjU1NDE2IC02LjI3MjQ3LC03LjM4MDM3IC0xMi40ODQzNiwtMTQuODE2OTIgLTE4LjY1MTIzLC0yMi4yOTUwNyAtNi4xNjY4OCwtNy40NzgxNSAtMTIuMjg4NzMsLTE0Ljk5NzkgLTE4LjM4MTExLC0yMi41NDQ2NSAtNi4wOTIzOCwtNy41NDY3NSAtMTIuMTU1MjgsLTE1LjEyMDUxIC0xOC4yMDQyNSwtMjIuNzA2NjcgLTYuMDQ4OTYsLTcuNTg2MTcgLTEyLjA4Mzk5LC0xNS4xODQ3NiAtMTguMTIwNjMsLTIyLjc4MTE2IC02LjAzNjY0LC03LjU5NjQgLTEyLjA3NDg5LC0xNS4xOTA2MiAtMTguMTMwMywtMjIuNzY4MDcgLTYuMDU1NDEsLTcuNTc3NDUgLTEyLjEyNzk3LC0xNS4xMzgxMyAtMTguMjMzMjMsLTIyLjY2NzQ0IC02LjEwNTI2LC03LjUyOTMyIC0xMi4yNDMyMiwtMTUuMDI3MjcgLTE4LjQyOTQyLC0yMi40NzkyNiAtNi4xODYxOSwtNy40NTE5OSAtMTIuNDIwNjMsLTE0Ljg1ODAzIC0xOC43MTg4NiwtMjIuMjAzNTIgLTAuMjI3OTEsMTUuMTY0OTYgLTAuMzY2NzQsMzAuMzMwOCAtMC40NDg5LDQ1LjQ5NzE4IC0wLjA4MjIsMTUuMTY2MzcgLTAuMTA3NjYsMzAuMzMzMjkgLTAuMTA4ODksNDUuNTAwNCAtMC4wMDEsMTUuMTY3MTIgMC4wMjE4LDMwLjMzNDQzIDAuMDM2Nyw0NS41MDE2MiAwLjAxNDksMTUuMTY3MTggMC4wMjE2LDMwLjMzNDIyIC0wLjAxMjIsNDUuNTAwOCAtOS4zMzMzMywwIC0xOC42NjY2NywwIC0yOCwwIC05LjMzMzMzLDAgLTE4LjY2NjY2LDAgLTI4LDAgMCwtMTEuOTE2NjcgMCwtMjMuODMzMzQgMCwtMzUuNzUgMCwtMTEuOTE2NjcgMCwtMjMuODMzMzQgMCwtMzUuNzUgMCwtMTEuOTE2NjYgMCwtMjMuODMzMzMgMCwtMzUuNzUgMCwtMTEuOTE2NjcgMCwtMzUuNzUgLTEwZS02LC0zNS43NSB6IFwiO1xudmFyIFRjaGFyID0gXCJjIDAsLTkuOTE2NjcgMCwtMTkuODMzMzQgMCwtMjkuNzUgMCwtOS45MTY2NyAwLC0xOS44MzMzNCAwLC0yOS43NSAwLC05LjkxNjY2IDAsLTE5LjgzMzMzIDAsLTI5Ljc1IDAsLTkuOTE2NjYgLTEuNDczOTQsLTE5LjgzMzMzIC0xLjQ3Mzk0LC0yOS43NSAwLC0xNS4zMzMzNCAtMjkuMTkyNzMsMCAtNDQuNTI2MDYsMCAtMTUuMzMzMzMsMCAtMzAuNjY2NjcsMCAtNDYsMCAwLC04IDAsLTE2IDAsLTI0LjAwMDAwMiAwLC04IDAsLTE2LjAwMDAwMSAwLC0yNC4wMDAwMDEgOS45MTY2NiwwIDE5LjgzMzMzLDAgMjkuNzUsMCA5LjkxNjY3LDAgMTkuODMzMzMsMCAyOS43NSwwIDkuOTE2NjYsMCAxOS44MzMzMywwIDI5Ljc1LDAgOS45MTY2NywwIDE5LjgzMzMzLDAgMjkuNzUsMCA5LjkxNjY2LDAgMTkuODMzMzMsMCAyOS43NSwwIDkuOTE2NjcsMCAxOS44MzMzMywwIDI5Ljc1LDAgOS45MTY2NiwwIDE5LjgzMzMzLDAgMjkuNzUsMCA5LjkxNjY3LDAgMTkuODMzMzMsMCAyOS43NSwwIDAsOCAwLDE2LjAwMDAwMSAwLDI0LjAwMDAwMSAwLDguMDAwMDAyIDAsMTYuMDAwMDAyIDAsMjQuMDAwMDAyIC0xNSwwIC0zMCwwIC00NSwwIC0xNSwwIC00NC44ODgwOSwtMTMuNTI1NjUgLTQ1LDEuNDczOTQgLTAuMDczNSw5Ljg1NTg4IC0wLjEwNjEzLDE4LjIzOTkgLTAuMTEyNDksMjguMDk5MjMgLTAuMDA2LDkuODU5MzIgMC4wMTM1LDE5LjcyMDAxIDAuMDQ1LDI5LjU4MTMzIDAuMDMxNSw5Ljg2MTMyIDAuMDc0NSwxOS43MjMyOSAwLjExNDMsMjkuNTg1MTcgMC4wMzk5LDkuODYxODggMC4wNzY1LDE5LjcyMzY3IDAuMDk1NCwyOS41ODQ2NyAwLjAxODksOS44NjA5OSAwLjAxOTksMTkuNzIxMTggLTAuMDExNywyOS41Nzk4NCAtMC4wMzE1LDkuODU4NjYgLTAuMDk1NiwxOS43MTU4IC0wLjIwNjg4LDI5LjU3MDY5IC0wLjExMTMxLDkuODU0ODggLTAuMjY5ODUsMTkuNzA3NTIgLTAuNDkwMzMsMjkuNTU3MTkgLTAuMjIwNDcsOS44NDk2NyAtMC41MDI4OSwxOS42OTYzNiAtMC44NjE5MywyOS41MzkzNyAtNy4zOTYyNCwtMC45OTk2NCAtMTguOTU2NTgsMC45NjcyNiAtMjkuNzA5MTIsMS44Mjk3MSAtMTAuNzUyNTMsMC44NjI0NSAtMjQuNDg1NTYsMi4wMjYwOSAtMjQuODYyMzEsLTQuNzk2OTYgLTAuNjYyMjMsLTExLjk5MzE0IC0wLjY2MjIzLC0yMS41NDM0OSAtMC40OTY2NywtMzAuNDgzMTQgMC4xNjU1NiwtOC45Mzk2NSAwLjQ5NjY3LC0xNy4yNjg2IDAuNDk2NjcsLTI2LjgxODk1IDAsLTkuNTUwMzUgMCwtMTkuMTAwNyAwLC0yOC42NTEwNSAwLC05LjU1MDM1IDAsLTE5LjEwMDY5IDAsLTI4LjY1MTA0IHpcIjtcbnZhciBTY2hhciA9IFwiYyAtMTYuNjU4NywtMi4zNDM4NyAtMzMuMzU5OTUsLTYuMjM1MTUgLTQ5LjI1MTM3LC0xMi4wODgzMSAtMTUuODkxNDIsLTUuODUzMTYgLTMwLjk3MzA1LC0xMy42NjgyIC00NC4zOTI2LC0yMy44NTk1NyAzLjczNzE1LC04LjAyNjM5IDcuOTYwNjUsLTE1Ljc5MzcxIDEyLjE5MTAxLC0yMy41NTIyOCA0LjIzMDM3LC03Ljc1ODU4IDguNDY3NTksLTE1LjUwODQyIDEyLjIzMjE5LC0yMy40OTk4NCA4LjgyOTAyLDYuNjkxNjIgMTguNDgyMzUsMTIuNzU2NjMgMjguNjc0MTIsMTcuODYzNDQgMTAuMTkxNzgsNS4xMDY4MSAyMC45MjIwNSw5LjI1NTQyIDMxLjkwNDg1LDEyLjExNDIxIDEwLjk4MjksMi44NTg4IDIyLjIxODMsNC40Mjc4IDMzLjQyMDYsNC4zNzUzOSAxMS4yMDIyLC0wLjA1MjQgMjIuMzcxMywtMS43MjYyMiAzMy4yMjEyLC01LjM1MzA0IDE0LjI4MDQsLTYuOTIxNjkgMTcuMzAzMywtMTkuNjQ0MzYgMTMuNDk0NiwtMzEuMTUwMjggLTMuODA4NywtMTEuNTA1OTEgLTE0LjQ0ODksLTIxLjc5NTA3IC0yNy40OTQ2LC0yMy44NDk3MiAtMTAuMTYyLC00LjA3Njg2IC0yMS4xMDUyLC03LjE0ODYzIC0zMi4xOTc1LC0xMC4xMTA5NyAtMTEuMDkyNCwtMi45NjIzNSAtMjIuMzM0LC01LjgxNTI2IC0zMy4wOTMxLC05LjQ1NDM5IC0xMC43NTkxLC0zLjYzOTE0IC0yMS4wMzU2MiwtOC4wNjQ0OSAtMzAuMTk3NzksLTE0LjE3MTcxIC05LjE2MjE2LC02LjEwNzIzIC0xNy4yMDk5NiwtMTMuODk2MzIgLTIzLjUxMTU4LC0yNC4yNjI5MyAtNS4xNzk0LC0xMC45NzA1OCAtNy40ODU0NCwtMjIuNzc1ODMgLTcuMzA1MTIsLTM0LjUyOTY2IDAuMTgwMzEsLTExLjc1MzgyIDIuODQ2OTgsLTIzLjQ1NjIgNy42MTI5OSwtMzQuMjIxMDMgNC43NjYwMSwtMTAuNzY0ODMgMTEuNjMxMzcsLTIwLjU5MjEgMjAuMjA5MDUsLTI4LjU5NTY5MiA4LjU3NzY5LC04LjAwMzU5MiAxOC44Njc3MSwtMTQuMTgzNTA2IDMwLjQ4MzA1LC0xNy42NTM2MjEgMTIuMTgyMSwtNC4wMjQxNDUgMjQuODc3NywtNi4zNTcwMDkgMzcuNjk4OCwtNy4xMjE5MDEgMTIuODIxMSwtMC43NjQ4OTIgMjUuNzY3NywwLjAzODE5IDM4LjQ1MTcsMi4yODU5MzMgMTIuNjg0MSwyLjI0Nzc0NCAyNS4xMDU1LDUuOTQwMTUzIDM2Ljg3NjUsMTAuOTUzOTE4IDExLjc3MDksNS4wMTM3NjQgMjIuODkxMiwxMS4zNDg4ODUgMzIuOTczLDE4Ljg4MjA1MyAtNC40NzQ1LDcuMDAzNjggLTguMjU4MSwxNC43MDcwOCAtMTIuMTg4NiwyMi4yMzQ5OSAtMy45MzA2LDcuNTI3OSAtOC4wMDgsMTQuODgwMzEgLTEzLjA3MDEsMjEuMTgxOTggLTE1LjA3MjgsLTExLjE1MTU4IC0zMy41MDA5LC0yMC41NDQ5MSAtNTIuNjExNCwtMjQuOTE3MjYgLTE5LjExMDUsLTQuMzcyMzUgLTM4LjkwMzQsLTMuNzIzNzIgLTU2LjcwNTgsNS4yMDg2MSAtMTAuMDMyMSw3LjAyNzg5IC0xMy4zMzQ4LDE4Ljg0Njk1IC0xMS4xNTYxLDI5LjUzNTk2IDIuMTc4NywxMC42ODkwMiA5LjgzODcsMjAuMjQ3OTkgMjEuNzMyLDIyLjc1NTcyIDEwLjEzMzUsNC40MzMwNyAyMS4wMjQxLDcuNjk0MDcgMzIuMTEwMSwxMC42OTAyMyAxMS4wODYxLDIuOTk2MTYgMjIuMzY3Niw1LjcyNzQ3IDMzLjI4Myw5LjEwMTE3IDEwLjkxNTUsMy4zNzM3IDIxLjQ2NDgsNy4zODk3OCAzMS4wODY1LDEyLjk1NTQ3IDkuNjIxNyw1LjU2NTY5IDE4LjMxNTcsMTIuNjgwOTkgMjUuNTIwNCwyMi4yNTMxMyA2LjYzMDEsOS43MTYzNSAxMC41NTg3LDIwLjg4MDIzIDEyLjAwNjQsMzIuNDEyMzEgMS40NDc2LDExLjUzMjA3IDAuNDE0MiwyMy40MzIzNCAtMi44Nzk3LDM0LjYyMTQ2IC0zLjI5MzksMTEuMTg5MTIgLTguODQ4NCwyMS42NjcwOSAtMTYuNDQzLDMwLjM1NDU3IC03LjU5NDYsOC42ODc0OSAtMTcuMjI5MywxNS41ODQ0OSAtMjguNjgzNywxOS42MTE2NiAtMTMuMDA5MSw1LjkyOTI4IC0yNy4wMTk3LDguNjM4NjkgLTQxLjI3MjgsOS42MzYwOCAtMTQuMjUzLDAuOTk3MzggLTI4Ljc0ODQsMC4yODI3NCAtNDIuNzI3MiwtMC42MzYwOCB6IFwiO1xuXG52YXIgc3ZnU3RyaW5nID0gXCJtIDAsMCBcIiArIEFjaGFyICtcIm0gMTI2LC0zMSBcIiArIE5jaGFyICsgXCJtIDM3NiwyNCBcIisgVGNoYXIgK1wibSAyNTAsMTIwIFwiKyBTY2hhcjtcblxuIFxuZnVuY3Rpb24gc3ZnVG9Qb2ludHMoc3ZnU3RyaW5nKSB7XG4gICAgdmFyIHBvaW50cyA9IFtdO1xuICAgIHZhciBlZGdlcyA9IE9iamVjdC5jcmVhdGUobnVsbCk7XG5cbiAgICB2YXIgYmVnaW5pbmdQYXRoO1xuXG4gICAgdmFyIFggPSAwO1xuICAgIHZhciBZID0gMDtcbiAgICB2YXIgbmJQb2ludHMgPSAwO1xuICAgIHZhciBwcmV2UG9pbnQ7XG5cbiAgICB2YXIgY29tbWFuZHMgPSBwYXJzZShzdmdTdHJpbmcpXG4gICAgZm9yICh2YXIgaT0wOyBpPGNvbW1hbmRzLmxlbmd0aDsgaSsrKXtcbiAgICAgICAgdmFyIGNvbW1hbmQgPSBjb21tYW5kc1tpXTtcbiAgICAgICAgc3dpdGNoIChjb21tYW5kWzBdKSB7XG4gICAgICAgICAgICBjYXNlIFwibVwiOlxuICAgICAgICAgICAgICAgIFggKz0gY29tbWFuZFsxXTtcbiAgICAgICAgICAgICAgICBZICs9IGNvbW1hbmRbMl07XG4gICAgICAgICAgICAgICAgcHJldlBvaW50ID0gdW5kZWZpbmVkO1xuICAgICAgICAgICAgICAgIGJlZ2luaW5nUGF0aCA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgY2FzZSBcIk1cIjpcbiAgICAgICAgICAgICAgICBYID0gY29tbWFuZFsxXTtcbiAgICAgICAgICAgICAgICBZID0gY29tbWFuZFsyXTtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYmVnaW5pbmdQYXRoID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgYnJlYWs7ICBcbiAgICAgICAgICAgIGNhc2UgXCJjXCI6XG4gICAgICAgICAgICAgICAgWCArPSBjb21tYW5kWzVdO1xuICAgICAgICAgICAgICAgIFkgKz0gY29tbWFuZFs2XTtcbiAgICAgICAgICAgICAgICBwb2ludHMucHVzaCh7aWQ6bmJQb2ludHMsIHg6WCwgeTpZfSk7XG5cbiAgICAgICAgICAgICAgICBpZiAocHJldlBvaW50ICE9IHVuZGVmaW5lZCkge1xuICAgICAgICAgICAgICAgICAgICBlZGdlc1twcmV2UG9pbnRdID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIG5iUG9pbnRzKys7XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBjYXNlIFwielwiOlxuICAgICAgICAgICAgICAgIGVkZ2VzW3ByZXZQb2ludF0gPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICBiZWdpbmluZ1BhdGggPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgcHJldlBvaW50ID0gdW5kZWZpbmVkO1xuICAgICAgICAgICAgICAgIGJyZWFrOyAgICBcbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4ge3BvaW50cyA6IHBvaW50cywgZWRnZXMgOiBlZGdlc307XG59XG5cbi8vIGluaXRpYWxpemUgcG9pbnRzXG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24obmJTdGFydFBvaW50cywgbmJSYW5kb21Qb2ludHMpe1xuICAgIHZhciBwb2ludHMgPSBbXTtcbiAgICB2YXIgZm9yY2VkRWRnZXM7XG4gICAgdmFyIGNpdHlTZXQ7XG5cbiAgICBpZiAodGV4dE1lc2gpe1xuXG4gICAgICAgIHZhciBteVRleHQgPSBzdmdUb1BvaW50cyhzdmdTdHJpbmcpO1xuICAgICAgICBwb2ludHMgPSBteVRleHQucG9pbnRzO1xuICAgICAgICBmb3JjZWRFZGdlcyA9IG15VGV4dC5lZGdlcztcbiAgICAgICAgY2l0eVNldCA9IHJhbmdlKDAsIHBvaW50cy5sZW5ndGgpO1xuXG4gICAgICAgIHZhciBzY2FsZVggPSAwLjU7XG4gICAgICAgIHZhciBzY2FsZVkgPSAwLjI7XG4gICAgICAgIHZhciBkZWx0YVggPSAwLjI1O1xuICAgICAgICB2YXIgZGVsdGFZID0gMC4zO1xuXG4gICAgICAgIC8vIHNjYWxlIHBvaW50cyB0byBbMCwxXSArIGRlbHRhXG4gICAgICAgIHZhciBtYXhYID0gTWF0aC5tYXguYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC54fSkpO1xuICAgICAgICB2YXIgbWluWCA9IE1hdGgubWluLmFwcGx5KE1hdGgsIHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIHAueH0pKTtcbiAgICAgICAgdmFyIG1heFkgPSBNYXRoLm1heC5hcHBseShNYXRoLCBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiBwLnl9KSk7XG4gICAgICAgIHZhciBtaW5ZID0gTWF0aC5taW4uYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC55fSkpO1xuICAgICAgICBwb2ludHMgPSBwb2ludHMubWFwKGZ1bmN0aW9uKHApe1xuICAgICAgICAgICAgdmFyIHggPSBzY2FsZVggKiAocC54LW1pblgpLyhtYXhYLW1pblgpICsgZGVsdGFYO1xuICAgICAgICAgICAgdmFyIHkgPSBzY2FsZVkgKiAocC55LW1pblkpLyhtYXhZLW1pblkpICsgZGVsdGFZO1xuICAgICAgICAgICAgdmFyIG5ld1BvaW50ID0gbmV3IFBvaW50KHgsIHkpO1xuICAgICAgICAgICAgbmV3UG9pbnQuaWQgPSBwLmlkO1xuXG4gICAgICAgICAgICByZXR1cm4gbmV3UG9pbnQ7XG4gICAgICAgIH0pO1xuXG4gICAgICAgIC8vIG9ubHkgYWRkIHJhbmRvbSBwb2ludHNcbiAgICAgICAgdmFyIG5iUG9pbnRzID0gcG9pbnRzLmxlbmd0aDtcbiAgICAgICAgZm9yKHZhciBpPTA7IGk8bmJSYW5kb21Qb2ludHM7ICsraSkge1xuXG4gICAgICAgICAgICB2YXIgeCA9IHJhbmRvbSgpO1xuICAgICAgICAgICAgdmFyIHkgPSByYW5kb20oKTtcblxuICAgICAgICAgICAgdmFyIG5ld1BvaW50ID0gbmV3IFBvaW50KHgsIHkpO1xuICAgICAgICAgICAgbmV3UG9pbnQuaWQgPSBuYlBvaW50cztcblxuICAgICAgICAgICAgcG9pbnRzLnB1c2gobmV3UG9pbnQpO1xuXG4gICAgICAgICAgICBuYlBvaW50cysrO1xuICAgICAgICB9XG5cbiAgICB9IGVsc2Uge1xuICAgICAgICAvL2FkZCByYW5kb20gcG9pbnRzXG5cbiAgICAgICAgdmFyIG5iUG9pbnRzID0gMDtcbiAgICAgICAgZm9yKHZhciBpPTA7IGk8bmJSYW5kb21Qb2ludHM7ICsraSkge1xuXG4gICAgICAgICAgICB2YXIgeCA9IHJhbmRvbSgpO1xuICAgICAgICAgICAgdmFyIHkgPSByYW5kb20oKTtcblxuICAgICAgICAgICAgdmFyIG5ld1BvaW50ID0gbmV3IFBvaW50KHgsIHkpO1xuICAgICAgICAgICAgbmV3UG9pbnQuaWQgPSBuYlBvaW50cztcblxuICAgICAgICAgICAgcG9pbnRzLnB1c2gobmV3UG9pbnQpO1xuICAgICAgICAgICAgXG4gICAgICAgICAgICBuYlBvaW50cysrO1xuICAgICAgICB9XG5cbiAgICAgICAgY2l0eVNldCA9IHJhbmdlKDAsIG5iQ2l0eSk7XG4gICAgICAgIGNvbnNvbGUubG9nKGNpdHlTZXQpO1xuICAgIH1cblxuXG4gICAgLy8gaW5pdGlhbGl6ZSBzdGFydCBwb2ludHNcbiAgICB2YXIgcG9zc2libGVTdGFydFBvaW50c0lkID0gW107XG5cbiAgICBpZiAobmJSYW5kb21Qb2ludHMgPT09IDApIHtcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBuYlN0YXJ0UG9pbnRzOyBpKyspe1xuICAgICAgICAgICAgcG9zc2libGVTdGFydFBvaW50c0lkLnB1c2goTWF0aC5mbG9vcihuYlJhbmRvbVBvaW50cyAqIHJhbmRvbSgpKSk7XG4gICAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IG5iU3RhcnRQb2ludHM7IGkrKyl7XG4gICAgICAgICAgICBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQucHVzaChNYXRoLmZsb29yKG5iUmFuZG9tUG9pbnRzICogcmFuZG9tKCkpKTtcbiAgICAgICAgfVxuICAgIH1cbiAgICBcblxuICAgIHJldHVybiB7XG4gICAgICAgIHRleHRNZXNoOiB0ZXh0TWVzaCxcbiAgICAgICAgcG9pbnRzOiBwb2ludHMsXG4gICAgICAgIGNpdHlTZXQ6IGNpdHlTZXQsXG4gICAgICAgIHBvc3NpYmxlU3RhcnRQb2ludHNJZDogcG9zc2libGVTdGFydFBvaW50c0lkLFxuICAgICAgICBuYlJhbmRvbVBvaW50czogbmJSYW5kb21Qb2ludHMsXG4gICAgICAgIGZvcmNlZEVkZ2VzOiBmb3JjZWRFZGdlc1xuICAgIH07XG59XG4iLCIndXNlIHN0cmljdCdcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbiAoY29udGFpbmVyKXtcblxuXHR2YXIgbW91c2UgPSB7XG5cdCAgICB4OiAwLFxuXHQgICAgeTogMFxuXHR9O1xuXG5cdGNvbnRhaW5lci5hZGRFdmVudExpc3RlbmVyKCAnbW91c2Vtb3ZlJywgZnVuY3Rpb24oZSl7XG5cdCAgICB2YXIgcmVjdCA9IGNvbnRhaW5lci5nZXRCb3VuZGluZ0NsaWVudFJlY3QoKTtcblxuXHQgICAgbW91c2UueCA9IChlLmNsaWVudFggLSByZWN0LmxlZnQgKSAvIHJlY3Qud2lkdGg7XG5cdCAgICBtb3VzZS55ID0gKGUuY2xpZW50WSAtIHJlY3QudG9wICkvIHJlY3QuaGVpZ2h0O1xuXHR9KTtcblxuXHRyZXR1cm4gbW91c2U7XG5cbn07XG4iLCIndXNlIHN0cmljdCdcblxuZnVuY3Rpb24gUG9pbnQoeCwgeSkge1xuICAgIHRoaXMuaWQgPSB1bmRlZmluZWQ7ICAgICAgICAgICAgICAgIFxuICAgIHRoaXMueCA9IHg7XG4gICAgdGhpcy55ID0geTtcbiAgICB0aGlzLm5leHRzID0gW107XG59XG5cbm1vZHVsZS5leHBvcnRzID0gUG9pbnQ7IiwiJ3VzZSBzdHJpY3QnXG5cbnZhciBhbnRGdW5jdGlvbiA9IHJlcXVpcmUoJy4vYW50LmpzJyk7XG5cbnZhciByYW5kb20gPSBNYXRoLnJhbmRvbTtcblxudmFyIFJBTkRPTU1WVCA9IDAuMDAzO1xudmFyIEFOVFNJWkUgPSAwLjAwMjtcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihjb250YWluZXIsIGluaXRWYXIpe1xuXHRcblx0aWYoIWNvbnRhaW5lcilcblx0XHR0aHJvdyBuZXcgVHlwZUVycm9yKCdNaXNzaW5nIGNvbnRhaW5lcicpO1xuXG5cdHZhciBlZGdlcyA9IGluaXRWYXIuZWRnZXM7XG5cdHZhciBwb3B1bGF0aW9uID0gaW5pdFZhci5wb3B1bGF0aW9uO1xuXHR2YXIgcG9pbnRzSW5mb3MgPSBpbml0VmFyLnBvaW50c0luZm9zO1xuXHR2YXIgbmJBbnRzID0gcG9wdWxhdGlvbi5sZW5ndGg7XG5cblx0dmFyIGFuaW1JRDtcblxuXHR2YXIgY2FudmFzTGlzdCA9IGRvY3VtZW50LmdldEVsZW1lbnRzQnlUYWdOYW1lKFwiY2FudmFzXCIpO1xuXHRcblx0aWYgKGNhbnZhc0xpc3QubGVuZ3RoID09PSAwKXtcblx0XHR2YXIgY2FudmFzID0gZG9jdW1lbnQuY3JlYXRlRWxlbWVudChcImNhbnZhc1wiKTtcblx0XHR2YXIgcmVjdCA9IGNvbnRhaW5lci5nZXRCb3VuZGluZ0NsaWVudFJlY3QoKTtcblx0XHRjYW52YXMud2lkdGggPSByZWN0LndpZHRoO1xuXHRcdGNhbnZhcy5oZWlnaHQgPSByZWN0LmhlaWdodDtcblx0XHRjYW52YXMuc3R5bGUuYmFja2dyb3VuZENvbG9yID0gXCJyZ2JhKDI1MCwgMjUwLCAyNTAsIDApXCI7IFxuXHRcdGNvbnRhaW5lci5hcHBlbmRDaGlsZChjYW52YXMpO1xuXHR9XG5cdGVsc2V7XG5cdFx0dmFyIGNhbnZhcyA9IGNhbnZhc0xpc3RbMF07XG5cdFx0Y29uc29sZS5sb2coJ0NBTlZBUycpO1xuXHR9XG5cdFxuXHR2YXIgY29udGV4dCA9IGNhbnZhcy5nZXRDb250ZXh0KFwiMmRcIik7XG5cdGNvbnRleHQuY2xlYXJSZWN0ICggMCAsIDAgLCBjYW52YXMud2lkdGgsIGNhbnZhcy5oZWlnaHQgKTtcblx0XG5cblx0ZnVuY3Rpb24gdGljaygpIHtcblx0XHR2YXIgdyA9IGNhbnZhcy53aWR0aDtcblx0XHR2YXIgaCA9IGNhbnZhcy5oZWlnaHQ7XG5cdFx0dmFyIG1vdXNlID0gW2xhc3RNb3VzZU1vdmVFdmVudC5jbGllbnRYL3csIGxhc3RNb3VzZU1vdmVFdmVudC5jbGllbnRZL2hdO1xuXHRcdGNvbnRleHQuc2V0VHJhbnNmb3JtKHcsIDAsIDAsIGgsIDAsIDApO1xuXHRcdGNvbnRleHQuZmlsbFN0eWxlID0gXCJyZ2JhKDI1MCwgMjUwLCAyNTAsIDAuNClcIjtcblx0XHRjb250ZXh0LmZpbGxSZWN0KDAsMCx3LGgpO1xuXG5cdFx0Ly8gZWRnZXNcblx0XHQvLyBjb250ZXh0LnN0cm9rZVN0eWxlID0gXCIjMDAwXCI7XG5cdFx0Ly8gZm9yKHZhciBpPTA7IGkgPCBlZGdlcy5sZW5ndGg7ICsraSkge1xuXHRcdC8vICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IDAuMDAwMTtcblx0XHQvLyAgICAgdmFyIGVkZ2UgPSBlZGdlc1tpXTtcblx0XHQvLyAgICAgLy8gaWYgKGVkZ2UucGhlcm9tb24gIT0gMCl7XG5cdFx0Ly8gICAgIC8vICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IE1hdGgubWluKDAuMDAwMDEgKiBlZGdlLnBoZXJvbW9uLCAwLjAxKTtcblx0XHQvLyAgICAgLy8gfSBlbHNlIHtcblx0XHQvLyAgICAgLy8gICAgIGNvbnRleHQubGluZVdpZHRoID0gMC4wMDAwMTtcblx0XHQvLyAgICAgLy8gfVxuXHRcdC8vICAgICBjb250ZXh0LmJlZ2luUGF0aCgpO1xuXHRcdC8vICAgICBjb250ZXh0Lm1vdmVUbyhwb2ludHNbZWRnZS5wdDEuaWRdLngsIHBvaW50c1tlZGdlLnB0MS5pZF0ueSk7XG5cdFx0Ly8gICAgIGNvbnRleHQubGluZVRvKHBvaW50c1tlZGdlLnB0Mi5pZF0ueCwgcG9pbnRzW2VkZ2UucHQyLmlkXS55KTtcblx0XHQvLyAgICAgY29udGV4dC5zdHJva2UoKTtcblx0XHQvLyB9XG5cblx0XHQvLyAvLyB2ZXJ0aWNlc1xuXHRcdC8vIGZvcih2YXIgaT0wOyBpPHBvaW50cy5sZW5ndGg7ICsraSkge1xuXHRcdC8vICAgICBjb250ZXh0LmJlZ2luUGF0aCgpXG5cdFx0Ly8gICAgIHZhciBwb2ludCA9IHBvaW50c1tpXTtcblx0XHQvLyAgICAgaWYgKGNpdHlTZXQuaGFzKHBvaW50LmlkKSkge1xuXHRcdC8vICAgICAgICAgY29udGV4dC5maWxsU3R5bGUgPSBcIiMwMTAxREZcIjtcblx0XHQvLyAgICAgICAgIGNvbnRleHQuYXJjKHBvaW50LngsIHBvaW50LnksIDAuMDA2LCAwLCAyKk1hdGguUEkpO1xuXHRcdC8vICAgICB9XG5cdFx0Ly8gICAgIGVsc2Uge1xuXHRcdC8vICAgICAgICAgY29udGV4dC5maWxsU3R5bGUgPSBcIiMwMDBcIjtcblx0XHQvLyAgICAgICAgIGNvbnRleHQuYXJjKHBvaW50c1tpXS54LCBwb2ludHNbaV0ueSwgMC4wMDMsIDAsIDIqTWF0aC5QSSk7XG5cdFx0Ly8gICAgIH1cblx0XHQvLyAgICAgY29udGV4dC5jbG9zZVBhdGgoKTtcblx0XHQvLyAgICAgY29udGV4dC5maWxsKCk7XG5cdFx0Ly8gfVxuXG5cdFx0Ly8gbW92ZSBhbnRzXG5cdFx0cG9wdWxhdGlvbi5mb3JFYWNoKGZ1bmN0aW9uKGFudCl7XG5cdFx0XHRhbnQudHJhbnNpdCgpO1xuXHRcdH0pO1xuXG5cdFx0Ly8gZm9yIChpID0gMDsgaSA8IG5iQW50czsgaSsrKSB7XG5cdFx0Ly8gICAgIHBvcHVsYXRpb25baV0udHJhbnNpdCgpO1xuXHRcdC8vIH1cblxuXHRcdC8vIHBoZXJvbW9uIGV2YXBvcmF0aW9uXG5cdFx0ZWRnZXMuZm9yRWFjaChmdW5jdGlvbihlZGdlKXtcblx0XHRcdGlmKGVkZ2UucGhlcm9tb24gPiAwKXtcblx0XHRcdFx0ZWRnZS5waGVyb21vbiAtPSAwLjAwMDE7XG5cdFx0XHR9XG5cdFx0fSk7XG5cblx0XHQvLyBmb3IgKGkgPSAwOyBpIDwgZWRnZXMubGVuZ3RoOyBpKyspIHtcblx0XHRcdFxuXG5cdFx0Ly8gYW50c1xuXHRcdHBvcHVsYXRpb24uZm9yRWFjaChmdW5jdGlvbihhbnQpe1xuXHRcdFx0Y29udGV4dC5iZWdpblBhdGgoKVxuXHRcdFx0dmFyIHggPSBhbnQueCArIFJBTkRPTU1WVCpyYW5kb20oKTtcblx0XHRcdHZhciB5ID0gYW50LnkgKyBSQU5ET01NVlQqcmFuZG9tKCk7XG5cblx0XHRcdGNvbnRleHQuZmlsbFN0eWxlID0gXCJibGFja1wiXG5cdFx0XHRjb250ZXh0LmZpbGxSZWN0KHgsIHksIEFOVFNJWkUsIEFOVFNJWkUpO1xuXHRcdFx0Y29udGV4dC5jbG9zZVBhdGgoKTtcblx0XHRcdGNvbnRleHQuZmlsbCgpO1xuXHRcdH0pXG5cdFx0Ly8gZm9yKHZhciBpPTA7IGk8cG9wdWxhdGlvbi5sZW5ndGg7ICsraSkge1xuXHRcdFx0XG5cdFx0Ly8gfVxuXG5cdH07XG5cdFxuXHR2YXIgbGFzdE1vdXNlTW92ZUV2ZW50ID0ge1xuXHRcdGNsaWVudFg6IDAsXG5cdFx0Y2xpZW50WTogMFxuXHR9O1xuXHRcblx0Y29udGFpbmVyLmFkZEV2ZW50TGlzdGVuZXIoJ21vdXNlbW92ZScsIGZ1bmN0aW9uKGUpe1xuXHRcdGxhc3RNb3VzZU1vdmVFdmVudCA9IGU7XG5cdH0pO1xuXHRcblx0dmFyIHBhdXNlZCA9IGZhbHNlO1xuXHRcblx0ZnVuY3Rpb24gdG9nZ2xlUGxheVBhdXNlKCl7XG5cdFx0cGF1c2VkID0gIXBhdXNlZDtcblx0XHRpZighcGF1c2VkKVxuXHRcdFx0YW5pbWF0ZSgpO1xuXHR9XG5cblx0ZnVuY3Rpb24gcmVzZXQoKXtcblx0XHRwb3B1bGF0aW9uID0gW107XG5cdFx0ZWRnZXMgPSBbXTtcblx0XHRwb2ludHNJbmZvcyA9IFtdO1xuXG5cdFx0Y2FuY2VsQW5pbWF0aW9uRnJhbWUoYW5pbUlEKTtcblx0fVxuXHRcblx0Ly8gY29udGFpbmVyLmFkZEV2ZW50TGlzdGVuZXIoJ2NsaWNrJywgdG9nZ2xlUGxheVBhdXNlKTtcblxuXHRmdW5jdGlvbiBhbmltYXRlKCl7XG5cdFx0dGljaygpO1xuXHRcdFxuXHRcdGlmKCFwYXVzZWQpXG5cdFx0XHRhbmltSUQgPSByZXF1ZXN0QW5pbWF0aW9uRnJhbWUoYW5pbWF0ZSk7XG5cdH1cblx0YW5pbWF0ZSgpO1xuXG5cblx0ZnVuY3Rpb24gc2V0QW50Q291bnQob3B0cyl7XG5cblx0XHR2YXIgcHJldmlvdXNDb3VudCA9IHBvcHVsYXRpb24ubGVuZ3RoO1xuXG5cdFx0Y29uc29sZS5sb2coJ25iQW50cyA6Jywgb3B0cy5uYkFudHMpO1xuXHRcdGNvbnNvbGUubG9nKCdwcmV2aW91c0NvdW50IDonLCBwcmV2aW91c0NvdW50KTtcblxuXG5cdFx0aWYgKG9wdHMubmJBbnRzID4gcHJldmlvdXNDb3VudCl7XG5cdFx0XHRjb25zb2xlLmxvZygncG9pbnRzSW5mb3MgOicsIHBvaW50c0luZm9zKTtcblx0XHRcdHZhciBBbnQgPSBhbnRGdW5jdGlvbihjb250YWluZXIsIHBvaW50c0luZm9zLCBvcHRzKTtcblxuXHRcdFx0Zm9yICh2YXIgaSA9IDA7IGkgPCBvcHRzLm5iQW50cyAtIHByZXZpb3VzQ291bnQ7IGkrKykge1xuXHRcdFx0XHRjb25zb2xlLmxvZygndGVzdCcpO1xuXHRcdFx0XHR2YXIgbmV3QW50ID0gbmV3IEFudChBbnQuZ2VuZXJhdGVSYW5kU3RhcnRQb2ludCgpKTtcblx0XHRcdFx0bmV3QW50LnNldERpcmVjdGlvbigpO1xuXHRcdFx0XHRwb3B1bGF0aW9uLnB1c2gobmV3QW50KTtcblx0XHRcdH1cblx0XHR9XG5cdFx0ZWxzZXtcblx0XHRcdHBvcHVsYXRpb24gPSBwb3B1bGF0aW9uLnNsaWNlKDAsIG9wdHMubmJBbnRzKTtcblx0XHRcdG5iQW50cyA9IHBvcHVsYXRpb24ubGVuZ3RoO1xuXHRcdFx0Y29uc29sZS5sb2coJ05iIEFudHMgOicsIHBvcHVsYXRpb24ubGVuZ3RoKTtcblx0XHR9XG5cdFx0Y29uc29sZS5sb2coJ05vdXZlYXUgY29tcHRlIDonLCBwb3B1bGF0aW9uLmxlbmd0aCk7XG5cdH1cblxuXHRmdW5jdGlvbiBtb2RpZnlBbnRzKG9wdHMpe1xuXHRcdHNldEFudENvdW50KG9wdHMpO1xuXG5cdFx0cG9wdWxhdGlvbi5mb3JFYWNoKGZ1bmN0aW9uKGFudCl7XG5cdFx0XHRhbnQudmVsb2NpdHkgPSBvcHRzLnZlbG9jaXR5O1xuXHRcdFx0YW50LndlaWdodCA9IG9wdHMud2VpZ2h0O1xuXHRcdFx0YW50LnJlcFNpemUgPSBvcHRzLnJlcFNpemU7XG5cdFx0XHRhbnQucmVwU3BlZWQgPSBvcHRzLnJlcFNwZWVkO1xuXHRcdH0pO1xuXHR9XG5cdFxuXHRyZXR1cm4ge1xuXHRcdHRvZ2dsZVBsYXlQYXVzZTogdG9nZ2xlUGxheVBhdXNlLFxuXHRcdHJlc2V0OiByZXNldCxcblx0XHQvLyBzaG91bGQgYmUgYSBnZXR0ZXIvc2V0dGVyLCBidXQgSUU4XG5cdFx0Z2V0QW50Q291bnQ6IGZ1bmN0aW9uKCl7XG5cdFx0XHRyZXR1cm4gcG9wdWxhdGlvbi5sZW5ndGg7XG5cdFx0fSxcblx0XHRzZXRBbnRDb3VudDogc2V0QW50Q291bnQsXG5cdFx0bW9kaWZ5QW50czogbW9kaWZ5QW50c1xuXHR9XG59XG4iLCIndXNlIHN0cmljdCc7XG5cbnZhciBzcXJ0ID0gTWF0aC5zcXJ0O1xudmFyIHBvdyA9IE1hdGgucG93O1xuXG5mdW5jdGlvbiBzaWduKHgpIHtcblx0cmV0dXJuIHggPyB4IDwgMCA/IC0xIDogMSA6IDA7XG59XG5cbmZ1bmN0aW9uIHJhbmdlKHN0YXJ0LCBjb3VudCkge1xuICAgIHJldHVybiBBcnJheS5hcHBseSgwLCBBcnJheShjb3VudCkpLm1hcChmdW5jdGlvbiAoZWxlbWVudCwgaW5kZXgpIHtcbiAgICBcdHJldHVybiBpbmRleCArIHN0YXJ0XG4gICAgfSk7XG59XG5cbmZ1bmN0aW9uIGRpc3RhbmNlKGEsIGIpe1xuXHRyZXR1cm4gc3FydChwb3coYS54IC0gYi54LCAyKSArIHBvdyhhLnkgLSBiLnksIDIpKTtcbn1cblxuZnVuY3Rpb24gbm9ybSh2KXtcblx0cmV0dXJuIHNxcnQocG93KHYueCwgMikgKyBwb3codi55LCAyKSk7XG59XG5cbm1vZHVsZS5leHBvcnRzID0ge1xuXHRzaWduOiBzaWduLFxuXHRyYW5nZTogcmFuZ2UsXG5cdGRpc3RhbmNlOiBkaXN0YW5jZSxcblx0bm9ybTogbm9ybVxufSIsIid1c2Ugc3RyaWN0J1xuXG5mdW5jdGlvbiBWZWN0b3IoeCwgeSkge1xuICAgIHRoaXMueCA9IHg7ICAgICAgICAgICAgICAgIFxuICAgIHRoaXMueSA9IHk7XG59XG5cblZlY3Rvci5wcm90b3R5cGUubm9ybSA9IGZ1bmN0aW9uKCl7XG5cdHJldHVybiBNYXRoLnNxcnQodGhpcy54ICogdGhpcy54ICsgdGhpcy55ICogdGhpcy55KTtcbn1cblxuVmVjdG9yLnByb3RvdHlwZS5ub3JtYWxpemUgPSBmdW5jdGlvbigpe1xuXHR2YXIgbm9ybSA9IHRoaXMubm9ybSgpO1xuXHR0aGlzLnggPSB0aGlzLnggLyBub3JtO1xuXHR0aGlzLnkgPSB0aGlzLnkgLyBub3JtO1xufVxuXG5cblxubW9kdWxlLmV4cG9ydHMgPSBWZWN0b3I7Il19
