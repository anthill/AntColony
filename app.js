(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({"/Users/Romain/Documents/Programmation/ants/AntColony/index.js":[function(require,module,exports){
'use strict';

var initRendering = require('./src/rendering.js');
var initializePoints = require('./src/initializePoints.js');
var createEdges = require('./src/createEdges.js');
// var initAnts = require('./src/initializeAnts');

module.exports = function init(containerElement, options){

	var render, pointsInfos, edges, population, pointsMap;


	function _init(containerElement, options){
		pointsInfos = initializePoints(options.nbStart, options.nbRand);
		edges = createEdges(pointsInfos.points);
		// population = options.nbAnts;
		// population = initAnts(containerElement, pointsInfos, options);
		pointsMap = {
			pointsInfos: pointsInfos,
			edges: edges
			// population: population
		};
		render = initRendering(containerElement, pointsMap, options);
	}

	_init(containerElement, options);

	return {
		togglePlayPause: function(){ render.togglePlayPause() },
		changeOptions: function(opts){
			render.modifyAnts(opts);
		},
		reset: function(opts){
			render.reset();

				// reset elements
			_init(containerElement, opts);
		}
	};
};
},{"./src/createEdges.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/createEdges.js","./src/initializePoints.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/initializePoints.js","./src/rendering.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/rendering.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js":[function(require,module,exports){
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
},{"robust-orientation":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js","simplicial-complex":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js":[function(require,module,exports){
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
},{"two-product":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js","two-sum":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js":[function(require,module,exports){
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
},{"robust-scale":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js","robust-subtract":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js","robust-sum":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js","two-product":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js":[function(require,module,exports){
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


},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js":[function(require,module,exports){
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

},{"bit-twiddle":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js","union-find":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js":[function(require,module,exports){
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

},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/triangulate.js":[function(require,module,exports){
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
},{"incremental-convex-hull":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js","uniq":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/parse-svg-path/index.js":[function(require,module,exports){

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

},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/ant.js":[function(require,module,exports){
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

            // console.log('smells1: ', possibleEdges);

            possibleEdges.splice(possibleEdges.indexOf(this.edge),1);

            // flip a coin and either take the smelliest path or a random one
            if (random() < 0.75){
                var smells = possibleEdges.map(function(e){return e.pheromon;});
                var index = smells.indexOf(Math.max.apply(Math, smells));
                this.edge = possibleEdges[index];
            } 
            else{
                this.edge = possibleEdges[floor(random()*possibleEdges.length)];
            }
                

            // set the destination point, being edge.pt1 or edge.pt2
            this.destination = (this.origin == this.edge.pt1) ? this.edge.pt2 : this.edge.pt1;

            this.direction.x = this.destination.x - this.origin.x; 
            this.direction.y = this.destination.y - this.origin.y;

            this.direction.normalize();
        },

        move: function(){
            // console.log('move');
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
                // console.log('reached');
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


},{"./mouse.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/mouse.js","./utilities.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/utilities.js","./vector.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/vector.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/antsGroup.js":[function(require,module,exports){
'use strict'

module.exports = function (Ant) {

	var nbAntsPerStep = 100;

	function createGroup(population){
		for (var i = 0; i < nbAntsPerStep; i++) {
			var newAnt = new Ant(Ant.generateRandStartPoint());
			newAnt.setDirection();
			population.push(newAnt);
		}

// 		console.log('Created Ants Group: \
// (+ ' + nbAntsPerStep + ') => ' + population.length);

		return population;
	}

	function removeGroup(population, nbDead){
		population = population.slice(0, population.length - nbDead);

// 		console.log('Removed Ants Group: \
// (- ' + nbAntsPerStep + ') => ' + population.length);

		return population;

	}

	return {
		create: createGroup,
		remove: removeGroup
	};

}
	
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/createEdges.js":[function(require,module,exports){
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
},{"./edge.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/edge.js","./initializePoints.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/initializePoints.js","./utilities.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/utilities.js","delaunay-triangulate":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/delaunay-triangulate/triangulate.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/edge.js":[function(require,module,exports){
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
    };

    if (this.distance === 0) console.log('ZERO !');
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
},{"./vector.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/vector.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/initializePoints.js":[function(require,module,exports){
'use strict';

var parse = require('parse-svg-path');

var range = require('./utilities.js').range;

var Point = require('./point.js');
var svgPath = require('./svgPath.js');

var random = Math.random;

var nbCity = 2;

var textMesh = true;

// Frame definition
var xInit = 0, yInit = 0;
var w = 1,
    h = 1;

var svgString = svgPath;

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
                X = command[1];
                Y = command[2];
                points.push({id:nbPoints, x:X, y:Y});

                if (prevPoint != undefined) {
                    edges[prevPoint] = nbPoints;
                }
                prevPoint = nbPoints;
                nbPoints++;
                break; 
            case "l":
                X = command[1];
                Y = command[2];
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
        var scaleY = 0.4;
        var deltaX = 0.25;
        var deltaY = 0.25;

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

    for (var i = 0; i < nbStartPoints; i++){
        possibleStartPointsId.push(Math.floor(nbPoints * random()));
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

},{"./point.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/point.js","./svgPath.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/svgPath.js","./utilities.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/utilities.js","parse-svg-path":"/Users/Romain/Documents/Programmation/ants/AntColony/node_modules/parse-svg-path/index.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/mouse.js":[function(require,module,exports){
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

},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/point.js":[function(require,module,exports){
'use strict'

function Point(x, y) {
    this.id = undefined;                
    this.x = x;
    this.y = y;
    this.nexts = [];
}

module.exports = Point;
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/rendering.js":[function(require,module,exports){
'use strict'

var antFunction = require('./ant.js');
var antsGroup = require('./antsGroup');

var random = Math.random;

var RANDOMMVT = 0.003;
var ANTSIZE = 0.002;

module.exports = function(container, pointsMap, options){

	if(!container)
		throw new TypeError('Missing container');

	// Ants variables
	var edges = pointsMap.edges;
	var objPopulationInitial = options.nbAnts;
	var objPopulation = objPopulationInitial;
	var pointsInfos = pointsMap.pointsInfos;
	var population = [];
	var nbAntsPerStep = 100;
	
	var Ant = antFunction(container, pointsInfos, options);
	antsGroup = antsGroup(Ant);

	// Animation variables
	var animID;
	var deltaTime;
	var FPSCount;
	var lastUpdate = performance.now();
	var FPSMonitor = document.querySelector('#FPS');
	var dTMonitor = document.querySelector('#dT');
	var warningMonitor = document.querySelector('#warning');
	var refreshTime = 0;
	var maxDeltaTime = 40;
	var FPSOverLimitCount = 0;
	var FPSUnderLimitCount = 0;

	// warning message disappears after 4 s
	warningMonitor.addEventListener("transitionend", function(){
		window.setTimeout(function(){
			warningMonitor.className = "invisible";
		}, 4000);
	})

	// Canvas
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
	

	function checkAntNumber(antNumber){
		if (antNumber < objPopulation - 50){
			FPSMonitor.style.color = "green";
			population = antsGroup.create(population);
		}	
		else if (antNumber > objPopulation){
			population = antsGroup.remove(population, antNumber - objPopulation);
			FPSMonitor.style.color = "red";
			warningMonitor.className = "visible";
		}
		else{
			FPSMonitor.style.color = "white";
			// warningMonitor.className = "invisible";
		}
	}

	function displayFPS(dT){
		FPSCount = (1000/dT).toFixed(2);
		var t = dT.toFixed(2);
		FPSMonitor.textContent = 'FPS : ' + FPSCount;  
		dTMonitor.textContent = 'nbAnts : ' + population.length;
		// dTMonitor.innerText = 'dT : ' + t + 'ms';
	}

	function tick() {
		var now = performance.now();
		deltaTime = now - lastUpdate;
		lastUpdate = now;
		refreshTime += deltaTime/1000; // in seconds

		// console.log('nbAnts', population.length);

		checkAntNumber(population.length);

		// display FPS info every 0.3 s
		if (refreshTime > 0.3){
			displayFPS(deltaTime);
			refreshTime = 0; 
		}

		// remove ants when frame rate is too low
		if (FPSOverLimitCount === 10) {
			objPopulation = objPopulation * maxDeltaTime / deltaTime;
			FPSOverLimitCount = 0;
		}

		while (FPSUnderLimitCount > 50 && objPopulation < objPopulationInitial) {
			objPopulation += 10;
		}

		// check duration of over/under framerate limit periods
		if (deltaTime > maxDeltaTime){
			FPSOverLimitCount++;
			FPSUnderLimitCount = 0;
		}
		else {
			FPSOverLimitCount = 0;
			FPSUnderLimitCount++;
		}

		// draw in canvas
		var w = canvas.width;
		var h = canvas.height;
		var mouse = [lastMouseMoveEvent.clientX/w, lastMouseMoveEvent.clientY/h];
		context.setTransform(w, 0, 0, h, 0, 0);
		context.fillStyle = "rgba(250, 250, 250, 0.4)";
		context.fillRect(0,0,w,h);

		// // edges
		// context.strokeStyle = "#000";
		// for(var i=0; i < edges.length; ++i) {
		//     context.lineWidth = 0.0001;
		//     var edge = edges[i];
		//     if (edge.pheromon != 0){
		//         context.lineWidth = Math.min(0.00001 * edge.pheromon, 0.01);
		//     } else {
		//         context.lineWidth = 0.00001;
		//     }
		//     context.beginPath();
		//     context.moveTo(pointsInfos.points[edge.pt1.id].x, pointsInfos.points[edge.pt1.id].y);
		//     context.lineTo(pointsInfos.points[edge.pt2.id].x, pointsInfos.points[edge.pt2.id].y);
		//     context.stroke();
		// }

		// // vertices
		// for(var i=0; i<pointsInfos.points.length; ++i) {
		//     context.beginPath()
		//     var point = pointsInfos.points[i];
		//     if (pointsInfos.citySet.indexOf(point.id) != -1){
		//         context.fillStyle = "#0101DF";
		//         context.arc(point.x, point.y, 0.006, 0, 2*Math.PI);
		//     }
		//     else {
		//         context.fillStyle = "#000";
		//         context.arc(pointsInfos.points[i].x, pointsInfos.points[i].y, 0.003, 0, 2*Math.PI);
		//     }
		//     context.closePath();
		//     context.fill();
		// }

		// move ants
		population.forEach(function(ant){
			ant.transit();
		});

		// pheromon evaporation
		edges.forEach(function(edge){
			if(edge.pheromon > 0){
				edge.pheromon -= 0.0001;
			}
		});

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

	function modifyAnts(opts){
		objPopulation = opts.nbAnts;

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
		modifyAnts: modifyAnts
	}
}

},{"./ant.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/ant.js","./antsGroup":"/Users/Romain/Documents/Programmation/ants/AntColony/src/antsGroup.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/svgPath.js":[function(require,module,exports){
'use strict';

var point1 = "m 18.25,19.5 c 18.25,15.7677688 15.0000951,12.75 11,12.75 c 6.99990488,12.75 3.75,15.7677688 3.75,19.5 c 3.75,23.2322312 6.99990488,26.25 11,26.25 c 15.0000951,26.25 18.25,23.2322312 18.25,19.5 z m 4.25,19.5 c 4.25,16.0525294 7.26810863,13.25 11,13.25 c 14.7318914,13.25 17.75,16.0525294 17.75,19.5 c 17.75,22.9474706 14.7318914,25.75 11,25.75 c 7.26810863,25.75 4.25,22.9474706 4.25,19.5 z";
var point2 = "m 89.25,8.5 c 89.25,4.21605787 85.5528716,0.75 81,0.75 c 76.4471284,0.75 72.75,4.21605787 72.75,8.5 c 72.75,12.7839421 76.4471284,16.25 81,16.25 c 85.5528716,16.25 89.25,12.7839421 89.25,8.5 z m 73.25,8.5 c 73.25,4.49967088 76.7163156,1.25 81,1.25 c 85.2836844,1.25 88.75,4.49967088 88.75,8.5 c 88.75,12.5003291 85.2836844,15.75 81,15.75 c 76.7163156,15.75 73.25,12.5003291 73.25,8.5 z";
var point3 = "m 160.25,11 c 160.25,5.33348375 155.208168,0.75 149,0.75 c 142.791832,0.75 137.75,5.33348375 137.75,11 c 137.75,16.6665162 142.791832,21.25 149,21.25 c 155.208168,21.25 160.25,16.6665162 160.25,11 z m 138.25,11 c 138.25,5.62082125 143.057903,1.25 149,1.25 c 154.942097,1.25 159.75,5.62082125 159.75,11 c 159.75,16.3791787 154.942097,20.75 149,20.75 c 143.057903,20.75 138.25,16.3791787 138.25,11 z";
var point4 = "m 160.25,76.5 c 160.25,72.7677688 157.000095,69.75 153,69.75 c 148.999905,69.75 145.75,72.7677688 145.75,76.5 c 145.75,80.2322312 148.999905,83.25 153,83.25 c 157.000095,83.25 160.25,80.2322312 160.25,76.5 z m 146.25,76.5 c 146.25,73.0525294 149.268109,70.25 153,70.25 c 156.731891,70.25 159.75,73.0525294 159.75,76.5 c 159.75,79.9474706 156.731891,82.75 153,82.75 c 149.268109,82.75 146.25,79.9474706 146.25,76.5 z";
var point5 = "m 95.25,76 c 95.25,70.3362795 90.4344065,65.75 84.5,65.75 c 78.5655935,65.75 73.75,70.3362795 73.75,76 c 73.75,81.6637205 78.5655935,86.25 84.5,86.25 c 90.4344065,86.25 95.25,81.6637205 95.25,76 z m 74.25,76 c 74.25,70.6180255 78.8364268,66.25 84.5,66.25 c 90.1635732,66.25 94.75,70.6180255 94.75,76 c 94.75,81.3819745 90.1635732,85.75 84.5,85.75 c 78.8364268,85.75 74.25,81.3819745 74.25,76 z";
var point6 = "m 20.25,75 c 20.25,70.9919338 16.7764995,67.75 12.5,67.75 c 8.22350046,67.75 4.75,70.9919338 4.75,75 c 4.75,79.0080662 8.22350046,82.25 12.5,82.25 c 16.7764995,82.25 20.25,79.0080662 20.25,75 z m 5.25,75 c 5.25,71.2760797 8.49222829,68.25 12.5,68.25 c 16.5077717,68.25 19.75,71.2760797 19.75,75 c 19.75,78.7239203 16.5077717,81.75 12.5,81.75 c 8.49222829,81.75 5.25,78.7239203 5.25,75 z";
var point7 = "m 23.25,139 c 23.25,132.786797 18.2132034,127.75 12,127.75 c 5.78679656,127.75 0.75,132.786797 0.75,139 c 0.75,145.213203 5.78679656,150.25 12,150.25 c 18.2132034,150.25 23.25,145.213203 23.25,139 z m 1.25,139 c 1.25,133.062939 6.06293894,128.25 12,128.25 c 17.9370611,128.25 22.75,133.062939 22.75,139 c 22.75,144.937061 17.9370611,149.75 12,149.75 c 6.06293894,149.75 1.25,144.937061 1.25,139 z";
var point8 = "m 95.2429209,133.141624 c 95.5537927,127.209836 90.7562698,122.141736 84.5337404,121.815627 c 78.311211,121.489518 73.0102086,126.028377 72.6993368,131.960165 c 72.388465,137.891953 77.1859879,142.960053 83.4085173,143.286162 c 89.6310467,143.612271 94.932049,139.073412 95.2429209,133.141624 z m 73.1986516,131.986333 c 73.4947711,126.336036 78.555388,122.003001 84.5075724,122.314942 c 90.4597567,122.626883 95.0397256,127.465159 94.7436061,133.115456 c 94.4474866,138.765754 89.3868696,143.098788 83.4346853,142.786847 c 77.4825009,142.474907 72.902532,137.63663 73.1986516,131.986333 z";
var point9 = "m 167.728023,135.062479 c 168.038327,129.141533 163.482117,124.089864 157.551662,123.779062 c 151.621208,123.468261 146.561914,128.016002 146.25161,133.936948 c 145.941307,139.857894 150.497517,144.909562 156.427971,145.220364 c 162.358426,145.531166 167.41772,140.983425 167.728023,135.062479 z m 146.750925,133.963116 c 147.046767,128.318121 151.870617,123.982018 157.525494,124.278377 c 163.180372,124.574737 167.52455,129.391317 167.228709,135.036311 c 166.932867,140.681305 162.109017,145.017409 156.454139,144.72105 c 150.799262,144.42469 146.455083,139.60811 146.750925,133.963116 z";
var lettresANT = "m 25.6185414,31.9447266 l 25.6185414,31.1003906 l 25.388268,30.4095703 l 24.9277211,30.0257812 l 24.5439321,29.5652344 l 24.0833852,29.1814453 l 23.6995961,28.7208984 l 23.0087758,28.490625 l 22.1453678,28.490625 l 21.1512711,28.490625 l 20.475768,28.490625 l 19.8617055,28.7208984 l 19.4779164,29.1814453 l 18.7870961,29.3349609 l 17.5478579,29.3349609 l 16.669907,29.3349609 l 16.0024971,29.3349609 l 15.3349802,29.3349609 l 14.5654164,29.3349609 l 13.9513539,29.5652344 l 13.5675649,30.0257812 l 12.8767446,30.1792969 l 12.0324086,30.1792969 l 11.4183461,30.4095703 l 11.0345571,30.8701172 l 10.3437367,31.0236328 l 9.72967424,31.2539062 l 9.34588518,31.7144531 l 8.8853383,32.0982422 l 8.50154924,32.5587891 l 8.04100237,32.9425781 l 7.88748674,33.6333984 l 7.6572133,34.2474609 l 7.19666643,34.63125 l 7.0431508,35.3220703 l 7.0431508,36.1664062 l 7.19666643,36.7804687 l 7.6572133,37.1642578 l 8.04100237,37.6248047 l 8.50154924,38.0085937 l 8.8853383,38.4691406 l 9.34588518,38.8529297 l 9.72967424,39.3134766 l 10.3437367,39.4669922 l 11.4725961,39.4669922 l 12.3122961,39.4669922 l 12.8743818,39.4669922 l 13.6496873,39.4669922 l 14.5398528,39.4669922 l 15.4097524,39.4669922 l 16.1005727,39.3134766 l 16.4843617,38.8529297 l 17.0984242,38.6226562 l 17.7892446,38.4691406 l 18.1730336,38.0085937 l 18.6335805,37.6248047 l 19.0173696,37.1642578 l 19.6314321,36.9339844 l 20.3222524,36.7804687 l 20.7060414,36.3199219 l 21.3201039,36.0896484 l 22.0109242,35.9361328 l 22.3947133,35.4755859 l 23.0087758,35.2453125 l 23.6995961,35.0917969 l 23.9298696,34.4777344 l 24.0833852,33.7869141 l 24.5439321,33.403125 l 24.9277211,32.9425781 l 25.388268,32.5587891 l 25.6185414,31.9447266 z m 37.5927602,43.9189453 l 37.2089711,44.3794922 l 36.7484242,44.7632813 l 36.3646352,45.2238281 l 35.6738149,45.3773438 l 35.0597524,45.6076172 l 34.6759633,46.0681641 l 33.985143,46.2216797 l 33.395172,46.2216797 l 32.4874219,46.2216797 l 31.4521352,46.2216797 l 30.8380727,46.0681641 l 30.4542836,45.6076172 l 29.7634633,45.3773438 l 29.1494008,45.2238281 l 28.7656117,44.7632813 l 28.3050649,44.3794922 l 27.9212758,43.9189453 l 27.2304555,43.6886719 l 26.616393,43.5351563 l 26.2326039,43.0746094 l 25.7720571,43.0746094 l 25.388268,43.5351563 l 24.6974477,43.6886719 l 23.8531117,43.6886719 l 23.2390492,43.9189453 l 22.8552602,44.3794922 l 22.1644399,44.5330078 l 21.5503774,44.7632813 l 21.1665883,45.2238281 l 20.475768,45.3773438 l 19.6314321,45.3773438 l 19.0173696,45.6076172 l 18.6335805,46.0681641 l 17.9427602,46.2216797 l 17.1279087,46.2216797 l 16.2540883,46.2216797 l 15.6400258,46.4519531 l 15.2562367,46.9125 l 14.7956899,46.9125 l 14.4119008,46.4519531 l 13.7210805,46.2216797 l 12.7926837,46.2216797 l 12.01698,46.2216797 l 11.36469,46.2216797 l 10.48491,46.2216797 l 9.4994008,46.2216797 l 8.8853383,46.0681641 l 8.50154924,45.6076172 l 7.81072893,45.3773438 l 7.19666643,45.2238281 l 6.81287737,44.7632813 l 6.12205705,44.5330078 l 5.50799455,44.3794922 l 5.12420549,43.9189453 l 4.66365862,43.5351563 l 4.27986955,43.0746094 l 3.81932268,42.6908203 l 3.43553362,42.2302734 l 2.97498674,41.8464844 l 2.59119768,41.3859375 l 2.1306508,41.0021484 l 1.97713518,40.3880859 l 1.74686174,39.6972656 l 1.28631487,39.3134766 l 1.13279924,38.6994141 l 1.13279924,37.8550781 l 0.902525804,37.1642578 l 0.441978929,36.7804688 l 0.288463304,36.1664063 l 0.288463304,35.2521433 l 0.288463304,34.3970875 l 0.288463304,33.6333984 c 0.288463304,33.5008705 0.441978929,32.9425781 0.441978929,32.9425781 l 0.902525804,32.5587891 l 1.13279924,31.9447266 l 1.13279924,31.1003906 l 1.28631487,30.4095703 l 1.74686174,30.0257813 l 1.97713518,29.4117188 l 2.1306508,28.7208984 l 2.59119768,28.3371094 l 2.97498674,27.8765625 l 3.43553362,27.4927734 l 3.81932268,27.0322266 l 4.27986955,26.6484375 l 4.66365862,26.1878906 l 5.27772112,25.9576172 l 6.12205705,25.9576172 l 6.81287737,25.8041016 l 7.19666643,25.3435547 l 7.81072893,25.1132813 l 8.65506487,25.1132813 l 9.34588518,24.9597656 l 9.72967424,24.4992188 l 10.3437367,24.2689453 l 10.9370585,24.2689453 l 11.6695294,24.2689453 l 12.1244404,24.2689453 l 12.7890925,24.2689453 l 13.7210805,24.2689453 c 13.9566084,24.2689453 14.4119008,24.1154297 14.4119008,24.1154297 l 14.7956899,23.6548828 l 15.4097524,23.4246094 l 15.8914929,23.4246094 l 16.3780342,23.4246094 l 16.7165068,23.4246094 l 17.3264618,23.4246094 l 17.9427602,23.4246094 c 18.1782881,23.4246094 18.6335805,23.2710938 18.6335805,23.2710938 l 19.0173696,22.8105469 l 19.6314321,22.5802734 l 20.1258875,22.5802734 l 20.6556616,22.5802734 l 21.3026338,22.5802734 l 22.4916732,22.5802734 l 23.4347831,22.5802734 l 23.8531117,22.5802734 l 24.5439321,22.4267578 l 24.9277211,21.9662109 l 25.388268,21.5824219 l 25.7720571,21.121875 l 26.3861196,20.8916016 l 27.0769399,20.7380859 l 27.3072133,20.1240234 l 27.3072133,19.3201067 l 27.3072133,18.4353516 l 27.0769399,17.7445313 l 26.616393,17.3607422 l 26.616393,16.9001953 l 27.0769399,16.5164063 l 27.0769399,16.0558594 l 26.616393,15.6720703 l 26.4628774,15.0580078 l 26.2326039,14.3671875 l 25.7720571,13.9833984 l 25.388268,13.5228516 l 24.9277211,13.1390625 l 24.5439321,12.6785156 l 23.8531117,12.4482422 l 23.0087758,12.4482422 l 22.3947133,12.2947266 l 22.0109242,11.8341797 l 21.5503774,11.8341797 l 21.1665883,12.2947266 l 20.475768,12.4482422 l 19.5356151,12.4482422 l 18.6789018,12.4482422 l 18.1006461,12.4482422 l 17.4906911,12.4482422 l 16.8692027,12.4482422 l 16.2962649,12.4482422 l 15.7558552,12.4482422 l 15.1644088,12.4482422 l 14.5654164,12.4482422 l 13.9513539,12.6785156 l 13.5675649,13.1390625 l 13.107018,13.5228516 l 12.7232289,13.9833984 l 12.2626821,14.3671875 l 12.1091664,15.0580078 l 12.1091664,15.7020899 l 12.1091664,16.694579 l 12.1091664,17.5910156 l 11.878893,18.2050781 l 11.4183461,18.5888672 l 11.0345571,19.0494141 l 10.5740102,19.4332031 l 10.1902211,19.89375 l 9.4994008,20.0472656 l 8.6444784,20.0472656 l 7.81072893,20.0472656 l 7.19666643,19.89375 l 6.81287737,19.4332031 l 6.12205705,19.2029297 l 5.50799455,19.0494141 l 5.12420549,18.5888672 l 4.66365862,18.2050781 l 4.51014299,17.5910156 l 4.51014299,16.7466797 l 4.27986955,16.0558594 l 3.81932268,15.6720703 l 3.66580705,15.0580078 l 3.66580705,14.2136719 l 3.81932268,13.5228516 l 4.27986955,13.1390625 l 4.51014299,12.525 l 4.66365862,11.8341797 l 5.12420549,11.4503906 l 5.50799455,10.9898438 l 5.96854143,10.6060547 l 6.35233049,10.1455078 l 6.96639299,9.91523438 l 7.6572133,9.76171875 l 8.04100237,9.30117188 l 8.65506487,9.07089844 l 9.34588518,8.91738281 l 9.72967424,8.45683594 l 10.3437367,8.2265625 l 10.9890621,8.2265625 l 11.4361001,8.2265625 l 11.9428075,8.2265625 l 12.2681583,8.2265625 l 12.9549793,8.2265625 l 13.7210805,8.2265625 l 14.4119008,8.07304688 l 14.7956899,7.6125 l 15.2562367,7.6125 l 15.6400258,8.07304688 l 16.2540883,8.2265625 l 16.960917,8.2265625 l 17.4915199,8.2265625 l 18.029098,8.2265625 l 18.7932683,8.2265625 l 19.6314321,8.2265625 l 20.3222524,8.45683594 l 20.7060414,8.91738281 l 21.3201039,9.07089844 l 22.0440828,9.07089844 l 22.879456,9.07089844 l 23.8531117,9.07089844 c 24.0950474,9.07089844 24.5439321,9.30117188 24.5439321,9.30117188 l 24.9277211,9.76171875 l 25.5417836,9.91523438 l 26.2326039,10.1455078 l 26.616393,10.6060547 l 27.0769399,10.9898438 l 27.4607289,11.4503906 l 27.9212758,11.8341797 l 28.3050649,12.2947266 l 28.7656117,12.6785156 l 29.1494008,13.1390625 l 29.6099477,13.5228516 l 29.8402211,14.2136719 l 29.9937367,14.8277344 l 30.4542836,15.2115234 l 30.6845571,15.9023438 l 30.6845571,16.5032774 c 30.6845571,16.7271532 30.6845571,16.951029 30.6845571,17.1749048 l 30.6845571,17.8765741 l 30.6845571,18.5420549 l 30.6845571,19.1652008 c 30.6845571,19.7091147 30.6845571,20.1240234 30.6845571,20.1240234 l 30.8380727,20.7380859 l 31.2986196,21.121875 l 31.528893,21.8126953 l 31.528893,22.6060737 c 31.528893,22.862869 31.528893,23.1196643 31.528893,23.3764596 l 31.528893,24.110588 c 31.528893,24.3770742 31.528893,24.6555846 31.528893,24.943475 l 31.528893,25.8874136 l 31.528893,26.7741689 c 31.528893,27.0744041 31.528893,27.3746394 31.528893,27.6748747 l 31.528893,28.5871138 l 31.528893,29.2621252 l 31.528893,30.0192512 c 31.528893,30.2644768 31.528893,30.5128063 31.528893,30.7612527 l 31.528893,31.5916534 l 31.528893,32.2756429 l 31.528893,33.1500361 l 31.528893,33.8138594 l 31.528893,34.5533746 c 31.528893,35.0310778 31.528893,35.3220703 31.528893,35.3220703 l 31.6824086,35.9361328 l 32.1429555,36.3199219 l 32.5267446,36.7804688 l 32.9872914,37.1642578 l 33.3710805,37.6248047 l 33.985143,37.7783203 l 34.6759633,38.0085938 l 35.0597524,38.4691406 l 35.6738149,38.6226563 l 36.3646352,38.8529297 l 36.7484242,39.3134766 l 37.3624867,39.4669922 l 38.0533071,39.6972656 l 38.2835805,40.3880859 l 38.2835805,41.0212323 l 38.2835805,42.0753937 l 38.2835805,42.9210938 l 38.0533071,43.5351563 l 37.5927602,43.9189453 z m 87.4181857,40.5416016 l 87.2646701,41.2324219 l 87.0343966,41.8464844 l 86.5738497,42.2302734 l 86.1900607,42.6908203 l 85.7295138,43.0746094 l 85.3457247,43.5351563 l 84.6549044,43.6886719 l 84.0408419,43.9189453 l 83.6570529,44.3794922 l 82.9662326,44.5330078 l 82.3828191,44.5330078 l 81.7615463,44.5330078 l 81.1493993,44.5330078 l 80.3991122,44.5330078 l 79.7698767,44.5330078 l 79.2674188,44.5330078 l 78.5992642,44.5330078 l 77.9471246,44.5330078 l 77.2597342,44.5330078 l 76.6694177,44.5330078 l 76.158997,44.5330078 l 75.584606,44.5330078 l 75.0341926,44.5330078 l 74.507757,44.5330078 l 74.2712009,44.5330078 l 73.7096039,44.5330078 l 72.8342013,44.5330078 c 72.868641,44.5000413 72.2201388,44.3794922 72.2201388,44.3794922 l 71.8363497,43.9189453 l 71.3758029,43.5351563 l 71.2222872,42.9210938 l 70.9920138,42.2302734 l 70.5314669,41.8464844 l 70.3779513,41.2324219 l 70.5314669,40.5416016 l 70.9920138,40.1578125 l 71.3758029,39.6972656 l 71.8363497,39.3134766 l 72.2201388,38.8529297 l 72.6806857,38.4691406 l 73.0644747,38.0085938 l 73.6785372,37.7783203 l 74.3693576,37.6248047 l 74.7531466,37.1642578 l 75.2136935,36.7804688 l 75.5974826,36.3199219 l 76.0580294,35.9361328 l 76.2883029,35.3220703 l 76.2883029,34.4777344 l 76.4418185,33.7869141 l 76.9023654,33.403125 l 77.1326388,32.7890625 l 77.1326388,32.0586517 l 77.1326388,31.1347947 l 77.1326388,30.0680557 l 77.1326388,29.2630136 l 77.1326388,28.3870225 l 77.1326388,27.5049475 l 77.1326388,26.6028315 l 77.1326388,25.7716644 l 77.1326388,24.8237403 l 77.1326388,23.9217138 l 77.1326388,23.0904572 l 77.1326388,22.1313494 l 77.1326388,21.1642789 l 77.1326388,20.2382747 l 77.1326388,19.5429215 l 77.1326388,18.9216488 l 77.1326388,18.1794139 l 77.1326388,17.5910156 l 77.2861544,16.9001953 l 77.7467013,16.5164063 l 77.9769747,15.9023438 l 77.7467013,15.2115234 l 77.2861544,14.8277344 l 77.1326388,14.2136719 l 76.9023654,13.5228516 l 76.4418185,13.1390625 l 76.0580294,12.6785156 l 75.3672091,12.4482422 l 74.7531466,12.2947266 l 74.3693576,11.8341797 l 73.6785372,11.6039063 l 73.0644747,11.4503906 l 72.6806857,10.9898438 l 71.9898654,10.7595703 l 71.2135081,10.7595703 l 70.288488,10.7595703 l 69.7711782,10.7595703 l 69.0539946,10.7595703 l 68.2968184,10.7595703 l 67.5658566,10.7595703 l 66.9045912,10.7595703 l 66.0795138,10.7595703 c 66.1163893,10.810201 65.4654513,10.9898438 65.4654513,10.9898438 l 65.0816622,11.4503906 l 64.3908419,11.6039063 l 63.7767794,11.8341797 l 63.6232638,12.525 l 63.3929904,13.1390625 l 62.9324435,13.5228516 l 62.5486544,13.9833984 l 62.0881076,14.3671875 l 61.7043185,14.8277344 l 61.2437716,15.2115234 l 61.090256,15.9023438 l 60.8599826,16.5164063 l 60.3994357,16.9001953 l 60.0156466,17.3607422 l 59.5550997,17.7445313 l 59.4015841,18.4353516 l 59.1713107,19.0494141 l 58.7107638,19.4332031 l 58.5572482,20.1240234 l 58.3269747,20.7380859 l 57.8664279,21.121875 l 57.7129122,21.8126953 l 57.4826388,22.4267578 l 57.0220919,22.8105469 l 56.8685763,23.5013672 l 56.8685763,23.9897102 l 56.8685763,24.698931 l 56.8685763,25.6548179 l 56.8685763,26.3251197 l 56.8685763,27.0902587 l 56.8685763,27.8004637 l 56.8685763,28.5673828 l 57.0220919,29.1814453 l 57.4826388,29.5652344 l 57.7129122,30.2560547 l 57.7129122,31.1003906 l 57.4826388,31.7144531 l 57.0220919,32.0982422 l 56.8685763,32.7890625 l 57.0220919,33.403125 l 57.4826388,33.7869141 l 57.7129122,34.4777344 l 57.7129122,35.3220703 l 57.8664279,35.9361328 l 58.3269747,36.3199219 l 58.7107638,36.7804688 l 59.1713107,37.1642578 l 59.5550997,37.6248047 l 60.1691622,37.7783203 l 60.8599826,38.0085938 l 61.2437716,38.4691406 l 61.7043185,38.8529297 l 62.0881076,39.3134766 l 62.5486544,39.6972656 l 62.7789279,40.3880859 l 62.9324435,41.0021484 l 63.3929904,41.3859375 l 63.3929904,41.8464844 l 62.9324435,42.2302734 l 62.7789279,42.9210938 l 62.5486544,43.5351563 l 61.8578341,43.6886719 l 61.2437716,43.9189453 l 60.8599826,44.3794922 l 60.1691622,44.5330078 l 59.3248263,44.5330078 l 58.7107638,44.7632813 l 58.3269747,45.2238281 l 57.8664279,45.2238281 l 57.4826388,44.7632813 l 57.0220919,44.7632813 l 56.6383029,45.2238281 l 55.9474826,45.3773438 l 55.1031466,45.3773438 l 54.4890841,45.2238281 l 54.1052951,44.7632813 l 53.6447482,44.7632813 l 53.2609591,45.2238281 l 52.5701388,45.3773438 l 51.4974911,45.3773438 l 50.5906332,45.3773438 l 49.9853754,45.3773438 l 49.1927951,45.3773438 c 49.2268402,45.3443039 48.5787326,45.2238281 48.5787326,45.2238281 l 48.1949435,44.7632813 l 47.5041232,44.5330078 l 46.8900607,44.3794922 l 46.5062716,43.9189453 l 46.0457247,43.5351563 l 45.6619357,43.0746094 l 45.2013888,42.6908203 l 45.0478732,42.0767578 l 45.0478732,41.2324219 l 45.2013888,40.5416016 l 45.6619357,40.1578125 l 46.0457247,39.6972656 l 46.5062716,39.3134766 l 46.8900607,38.8529297 l 47.5041232,38.6226563 l 48.1949435,38.4691406 l 48.5787326,38.0085938 l 49.0392794,37.6248047 l 49.4230685,37.1642578 l 49.8836154,36.7804688 l 50.2674044,36.3199219 l 50.7279513,35.9361328 l 50.9582247,35.3220703 l 50.9582247,34.6998663 l 50.9582247,33.9150441 l 50.9582247,32.9400109 l 50.9582247,32.3075544 l 50.9582247,31.4683352 l 50.9582247,30.5891232 l 50.9582247,29.8138742 l 50.9582247,29.023684 l 50.9582247,28.1238047 l 50.9582247,27.0039211 l 50.9582247,26.034375 c 50.9582249,25.8460125 51.1117404,25.3435547 51.1117404,25.3435547 l 51.5722872,24.9597656 l 51.8025607,24.3457031 l 51.8025607,23.5013672 l 51.5722872,22.8105469 l 51.1117404,22.4267578 l 50.9582247,21.8126953 l 50.9582247,21.0290912 l 50.9582247,20.3358854 l 50.9582247,19.6347167 l 50.9582247,18.9575258 l 50.9582247,18.0005653 l 50.9582247,17.3553148 l 50.9582247,16.7466797 l 50.7279513,16.0558594 l 50.2674044,15.6720703 l 50.1138888,15.0580078 l 49.8836154,14.3671875 l 49.4230685,13.9833984 l 49.0392794,13.5228516 l 48.5787326,13.1390625 l 48.1949435,12.6785156 l 47.7343966,12.2947266 l 47.3506076,11.8341797 l 46.8900607,11.4503906 l 46.5062716,10.9898438 l 45.8154513,10.7595703 l 45.2013888,10.6060547 l 45.2013888,10.1455078 l 45.6619357,9.76171875 l 45.8922091,9.14765625 l 45.8922091,8.30332031 l 46.0457247,7.6125 l 46.5062716,7.22871094 l 46.8900607,6.76816406 l 47.3506076,6.384375 l 47.7343966,5.92382813 l 48.3484591,5.69355469 l 49.0392794,5.54003906 l 49.4230685,5.07949219 l 50.037131,4.84921875 l 50.7279513,5.07949219 l 51.1117404,5.54003906 l 51.7258029,5.69355469 l 52.5701388,5.69355469 l 53.2609591,5.92382813 l 53.6447482,6.384375 l 54.2588107,6.53789063 l 54.949631,6.76816406 l 55.3334201,7.22871094 l 55.7939669,7.6125 l 56.177756,8.07304688 l 56.6383029,8.45683594 l 57.0220919,8.91738281 l 57.4826388,9.30117188 l 57.8664279,9.76171875 l 58.3269747,10.1455078 l 58.7107638,10.6060547 l 59.1713107,10.6060547 l 59.5550997,10.1455078 l 60.1691622,9.91523438 l 61.0134982,9.91523438 l 61.7043185,9.76171875 l 62.0881076,9.30117188 l 62.5486544,8.91738281 l 62.9324435,8.45683594 l 63.3929904,8.07304688 l 63.7767794,7.6125 l 64.3908419,7.38222656 l 65.0816622,7.22871094 l 65.4654513,6.76816406 l 66.0795138,6.53789063 l 66.6456679,6.53789063 l 67.302102,6.53789063 l 68.0591887,6.53789063 l 68.6125216,6.53789063 l 69.3033419,6.384375 l 69.687131,5.92382813 l 70.3011935,5.69355469 l 70.9908197,5.69355469 l 71.5481222,5.69355469 l 72.0436015,5.69355469 l 72.7858364,5.69355469 l 73.6785372,5.69355469 l 74.3693576,5.92382813 l 74.7531466,6.384375 l 75.3672091,6.53789063 l 76.0580294,6.76816406 l 76.4418185,7.22871094 l 77.055881,7.38222656 l 77.7467013,7.6125 l 78.1304904,8.07304688 l 78.5910372,8.45683594 l 78.9748263,8.91738281 l 79.4353732,9.30117188 l 79.8191622,9.76171875 l 80.2797091,10.1455078 l 80.5099826,10.8363281 l 80.6634982,11.4503906 l 81.1240451,11.8341797 l 81.3543185,12.525 l 81.5078341,13.1390625 l 81.968381,13.5228516 l 82.1986544,14.2136719 l 82.1986544,14.9400989 l 82.1986544,15.9791025 l 82.1986544,16.7466797 l 81.968381,17.3607422 l 81.5078341,17.7445313 l 81.5078341,18.2050781 l 81.968381,18.5888672 l 81.968381,19.0494141 l 81.5078341,19.4332031 l 81.3543185,20.1240234 l 81.3543185,20.9653893 l 81.3543185,21.8126953 c 81.3543193,22.4267582 81.5078341,22.4267578 81.5078341,22.4267578 l 81.968381,22.8105469 l 82.1986544,23.5013672 l 82.1986544,24.3423076 l 82.1986544,25.2192828 l 82.1986544,25.9295772 l 82.1986544,27.0713807 l 82.1986544,27.6846908 l 82.1986544,28.5673828 c 82.1986529,29.1814445 82.3521701,29.1814453 82.3521701,29.1814453 l 82.8127169,29.5652344 l 83.0429904,30.2560547 l 83.0429904,31.1003906 l 82.8127169,31.7144531 l 82.3521701,32.0982422 l 82.1986544,32.7890625 l 82.3521701,33.403125 l 82.8127169,33.7869141 l 83.0429904,34.4777344 l 83.196506,35.0917969 l 83.6570529,35.4755859 l 84.0408419,35.9361328 l 84.6549044,36.0896484 l 85.3457247,36.3199219 l 85.7295138,36.7804688 l 86.3435763,36.9339844 l 87.0343966,37.1642578 l 87.2646701,37.8550781 l 87.4181857,38.4691406 l 87.8787326,38.8529297 l 88.109006,39.54375 l 87.8787326,40.1578125 l 87.4181857,40.5416016 z m 130.488924,32.0982422 l 130.335408,32.7890625 l 130.335408,33.6886606 l 130.335408,34.4777344 l 130.105135,35.0917969 l 129.644588,35.4755859 l 129.491072,36.1664063 l 129.491072,36.9742631 l 129.491072,37.8550781 l 129.260799,38.4691406 l 128.800252,38.8529297 l 128.646736,39.54375 l 128.416463,40.1578125 l 127.955916,40.5416016 l 127.572127,41.0021484 l 127.11158,41.3859375 l 126.727791,41.8464844 l 126.267244,42.2302734 l 125.883455,42.6908203 l 125.422908,43.0746094 l 125.039119,43.5351563 l 124.578572,43.9189453 l 124.194783,44.3794922 l 123.503963,44.5330078 l 122.8899,44.7632813 l 122.506111,45.2238281 l 121.815291,45.3773438 l 120.970955,45.3773438 l 120.356892,45.6076172 l 119.973103,46.0681641 l 119.282283,46.2216797 l 118.624185,46.2216797 l 117.941566,46.2216797 l 117.118507,46.2216797 l 116.442043,46.2216797 l 115.707548,46.2216797 l 115.060603,46.2216797 l 114.446541,46.0681641 l 114.062752,45.6076172 l 113.371932,45.3773438 l 112.527596,45.3773438 l 111.913533,45.2238281 l 111.529744,44.7632813 l 110.838924,44.5330078 l 110.224861,44.3794922 l 109.841072,43.9189453 l 109.380525,43.5351563 l 108.996736,43.0746094 l 108.536189,42.6908203 l 108.1524,42.2302734 l 107.691853,41.8464844 l 107.308064,41.3859375 l 106.847517,41.0021484 l 106.463728,40.5416016 l 106.003182,40.1578125 l 105.619392,39.6972656 l 105.158846,39.3134766 l 105.00533,38.6994141 l 104.775057,38.0085938 l 104.31451,37.6248047 l 104.160994,37.0107422 l 104.160994,36.1664063 l 103.930721,35.4755859 l 103.470174,35.0917969 l 103.316658,34.4777344 l 103.316658,33.7440534 l 103.316658,32.7890625 c 103.316658,32.4311707 103.470174,32.0982422 103.470174,32.0982422 l 103.930721,31.7144531 l 104.160994,31.1003906 l 104.160994,30.3962639 l 104.160994,29.5150941 l 104.160994,28.4306574 l 104.160994,27.7230469 c 104.160994,27.5442917 104.31451,27.0322266 104.31451,27.0322266 l 104.775057,26.6484375 l 105.00533,26.034375 l 104.775057,25.3435547 l 104.31451,24.9597656 l 104.160994,24.3457031 l 104.31451,23.6548828 l 104.775057,23.2710938 l 105.00533,22.6570313 l 105.00533,21.5411634 l 105.00533,20.6271416 l 105.00533,19.9922422 l 105.00533,19.2567085 l 105.00533,18.4955965 l 105.00533,17.6477584 l 105.00533,16.8000002 l 105.00533,15.9023438 l 104.775057,15.2115234 l 104.31451,14.8277344 l 104.160994,14.2136719 l 103.930721,13.5228516 l 103.470174,13.1390625 l 103.086385,12.6785156 l 102.395564,12.4482422 l 101.781502,12.2947266 l 101.397713,11.8341797 l 100.706892,11.6039063 l 100.053357,11.6039063 l 99.4288483,11.6039063 l 98.8200067,11.6039063 l 98.1738846,11.6039063 l 97.5598221,11.4503906 l 97.1760331,10.9898438 l 96.7154862,10.6060547 l 96.3316971,10.1455078 l 95.8711503,9.76171875 l 95.4873612,9.30117188 l 95.0268143,8.91738281 l 94.8732987,8.30332031 l 95.0268143,7.6125 l 95.4873612,7.22871094 l 95.8711503,6.76816406 l 96.3316971,6.384375 l 96.7154862,5.92382813 l 97.1760331,5.54003906 l 97.5598221,5.07949219 l 98.1738846,4.84921875 l 98.7337601,4.84921875 l 99.2183875,4.84921875 l 99.7019758,4.84921875 l 100.17853,4.84921875 l 100.604088,4.84921875 l 101.403807,4.84921875 l 102.002657,4.84921875 l 102.395564,4.84921875 l 103.086385,5.07949219 l 103.470174,5.54003906 l 103.930721,5.54003906 l 104.31451,5.07949219 l 104.775057,4.69570313 l 105.00533,4.08164063 l 105.158846,3.39082031 l 105.619392,3.00703125 l 105.849666,2.39296875 l 106.003182,1.70214844 l 106.617244,1.471875 l 107.308064,1.31835938 l 107.538338,0.704296875 l 107.691853,0.0134765625 l 108.1524,0.0134765625 l 108.382674,0.704296875 l 108.536189,1.31835938 l 108.996736,1.70214844 l 109.380525,2.16269531 l 109.841072,2.54648438 l 110.071346,3.23730469 l 110.071346,4.08164063 l 110.224861,4.69570313 l 110.685408,5.07949219 l 111.069197,5.54003906 l 111.68326,5.69355469 l 112.37408,5.92382813 l 112.757869,6.384375 l 113.218416,6.384375 l 113.602205,5.92382813 l 114.216267,5.69355469 l 115.004627,5.69355469 l 115.904939,5.69355469 l 116.59576,5.92382813 l 116.979549,6.384375 l 117.593611,6.53789063 l 118.251782,6.53789063 l 118.732093,6.53789063 l 119.389054,6.53789063 c 119.610092,6.53789063 119.83113,6.53789063 120.052169,6.53789063 c 120.203765,6.53789063 120.381803,6.53789063 120.572446,6.53789063 c 120.764922,6.53789063 120.957398,6.53789063 121.149874,6.53789063 c 121.352135,6.53789063 121.635686,6.53789063 121.902513,6.53789063 l 122.659627,6.53789063 c 122.89757,6.53789063 123.350447,6.76816406 123.350447,6.76816406 l 123.734236,7.22871094 l 124.194783,7.6125 l 124.578572,8.07304688 l 125.039119,8.45683594 l 125.269392,9.14765625 l 125.039119,9.76171875 l 124.578572,10.1455078 l 124.194783,10.6060547 l 123.734236,10.9898438 l 123.350447,11.4503906 l 122.659627,11.6039063 l 121.9882,11.6039063 l 121.349384,11.6039063 l 120.780988,11.6039063 c 120.593601,11.6039063 120.406214,11.6039063 120.218827,11.6039063 l 119.442448,11.6039063 l 118.774137,11.6039063 l 118.176246,11.6039063 l 117.520245,11.6039063 l 116.836667,11.6039063 l 116.238776,11.6039063 c 116.051863,11.6039063 115.868649,11.6039063 115.690922,11.6039063 c 115.490622,11.6039063 115.297292,11.6039063 115.113494,11.6039063 l 114.536066,11.6039063 c 114.374907,11.6039063 114.226441,11.6039063 114.093322,11.6039063 c 113.645967,11.6039063 113.371932,11.6039063 113.371932,11.6039063 l 112.757869,11.8341797 l 112.37408,12.2947266 l 111.913533,12.6785156 l 111.529744,13.1390625 l 111.069197,13.5228516 l 110.685408,13.9833984 l 110.224861,14.3671875 l 110.071346,15.0580078 l 110.071346,15.9902896 l 110.071346,16.7466797 c 110.071345,17.2189228 110.224861,17.3607422 110.224861,17.3607422 l 110.685408,17.7445313 l 110.685408,18.2050781 l 110.224861,18.5888672 l 110.071346,19.2796875 l 110.071346,20.1730483 c 110.071346,20.3950552 110.071346,20.6346766 110.071346,20.8877198 l 110.071346,21.6569848 l 110.071346,22.3905202 l 110.071346,23.1873618 l 110.071346,23.9709346 l 110.071346,24.6535534 l 110.071346,25.4737351 l 110.071346,26.265381 l 110.071346,26.9132294 l 110.071346,27.5442921 c 110.071346,28.1742542 110.071346,28.5673828 110.071346,28.5673828 l 110.224861,29.1814453 l 110.685408,29.5652344 l 110.915682,30.2560547 l 110.915682,31.1003906 l 111.069197,31.7144531 l 111.529744,32.0982422 l 111.760017,32.7890625 l 111.760017,33.6333984 l 111.913533,34.2474609 l 112.37408,34.63125 l 112.604353,35.3220703 l 112.757869,35.9361328 l 113.371932,36.0896484 l 114.062752,36.3199219 l 114.446541,36.7804688 l 115.060603,36.9339844 l 115.742318,36.9339844 l 116.226866,36.9339844 l 116.683917,36.9339844 l 117.041852,36.9339844 l 117.727348,36.9339844 l 118.437947,36.9339844 l 119.128767,36.7804688 l 119.512557,36.3199219 l 119.973103,35.9361328 l 120.356892,35.4755859 l 120.817439,35.0917969 l 121.047713,34.4777344 l 121.201228,33.7869141 l 121.661775,33.403125 l 121.892049,32.7890625 l 121.892049,31.9447266 l 122.045564,31.2539063 l 122.506111,30.8701172 l 122.736385,30.2560547 l 122.8899,29.5652344 l 123.350447,29.1814453 l 123.734236,28.7208984 l 124.348299,28.490625 l 125.039119,28.3371094 l 125.422908,27.8765625 l 126.036971,27.6462891 l 126.881307,27.6462891 l 127.572127,27.8765625 l 127.955916,28.3371094 l 128.569978,28.490625 l 129.260799,28.7208984 l 129.644588,29.1814453 l 130.105135,29.5652344 l 130.335408,30.2560547 l 130.488924,30.8701172 l 130.949471,31.2539063 l 130.949471,31.7144531 l 130.488924,32.0982422 z ";
var lettreS = "m 138.811517,35.1061011 c 138.63781,35.7365111 138.606564,35.8010614 138.198455,36.2091704 c 138.056753,36.3508722 137.999155,37.0215928 137.948798,37.2230199 c 137.869537,37.5400622 137.926411,37.8292331 137.948798,38.0775475 c 137.983078,38.4577806 137.927766,38.9415367 138.126037,39.2217479 c 138.240713,39.3838164 138.438596,39.433835 138.560544,39.6201964 c 138.747066,39.9052408 138.79591,40.4584527 138.951599,40.833933 c 139.061131,41.0980922 139.263918,41.2045958 139.429556,41.3702343 c 139.592672,41.5333495 139.626412,41.9345363 139.622169,42.4245876 c 139.61832,42.8690625 139.583226,43.386641 139.585076,43.8659947 c 139.586579,44.2554594 139.612471,44.6196912 139.699324,44.8989806 c 139.798518,45.2179557 140.058533,45.3368515 140.249382,45.5162438 c 140.468372,45.7220877 140.667698,45.9527405 140.883714,46.2064196 c 141.287987,46.6811791 141.363563,46.8042662 141.758733,46.8889452 c 142.153903,46.9736241 142.389571,46.8889452 142.5773,46.8889452 c 142.919018,46.8889452 143.062586,46.5447924 143.170055,46.4373229 c 143.5244,46.0829776 144.30261,46.0849327 144.456107,45.9576928 c 144.609603,45.8304528 144.538252,45.788898 144.680169,45.6469811 c 144.867391,45.4597592 145.3408,45.3492395 145.921047,45.2919329 c 146.416879,45.2429633 146.990723,45.23285 147.530673,45.2469365 c 148.234599,45.2653008 148.880917,45.3247952 149.221661,45.3929442 c 149.793527,45.5073172 149.592293,45.9406794 150.311344,46.0844896 c 150.590615,46.1403437 150.983508,46.1328596 151.385705,46.1339833 c 151.740461,46.1349745 152.102455,46.1426626 152.400101,46.2064196 c 153.035195,46.3424596 152.843259,46.6989649 153.303797,46.8553563 c 153.740109,47.0035208 154.126659,47.0607177 154.462767,47.0525307 c 154.848858,47.0431263 155.168389,46.9474451 155.420329,46.8042662 c 155.701153,46.6446723 155.856582,46.3018852 156.238444,46.2064196 c 156.63755,46.106643 157.062463,46.0773864 157.501193,46.0799003 c 157.962514,46.0825436 158.439111,46.1203133 158.917047,46.1481609 c 159.395158,46.1760187 159.874609,46.193947 160.341445,46.1568477 c 160.695624,46.1287011 161.042543,46.0688809 161.376107,45.9576928 c 161.556165,45.8976734 161.595501,45.6542159 161.855957,45.5162438 c 162.214357,45.3263872 162.76045,45.2045722 162.999127,45.1250133 c 163.301469,45.0242326 163.254211,44.8176843 163.542798,44.673391 c 163.811921,44.5388293 164.198751,44.5040294 164.537469,44.3628991 c 164.876187,44.2217687 165.245628,43.7958718 165.59596,43.4455407 c 165.694074,43.3474263 166.089427,42.951577 166.414031,42.5846337 c 166.738635,42.2176903 166.649398,41.6958381 166.82381,41.3702343 c 166.942587,41.1484928 167.231728,41.0717202 167.346164,40.833933 c 167.514797,40.4835286 167.470564,39.9584339 167.642542,39.6201964 c 167.752332,39.4042681 168.064756,39.3304646 168.164732,39.1262328 c 168.337576,38.7731432 168.382235,38.2192425 168.517395,37.9548375 c 168.663802,37.6684294 168.946353,37.5366951 169.011358,37.4066857 c 169.16199,37.1054214 169.196418,36.5836194 169.216001,36.0534301 c 169.23683,35.4895055 169.240866,34.9160924 169.350074,34.5884691 c 169.469774,34.2293693 169.81581,34.2906902 169.971058,33.7685027 c 170.126306,33.2463152 169.678613,32.7676354 169.40929,32.4983129 c 169.209305,32.2983273 169.341035,31.4437653 168.57385,30.8470659 c 168.319811,30.6494805 168.456555,30.0086667 168.164732,29.6053234 c 167.993437,29.3685691 167.724169,29.2856739 167.642542,29.0407946 c 167.455896,28.4808568 167.609275,28.4279622 167.346164,27.982304 c 167.124194,27.6063293 166.62278,27.1459759 166.287673,26.8554088 c 165.944635,26.5579643 165.711054,26.2484947 165.454331,26.0064531 c 165.144628,25.7144613 165.045061,25.484594 164.678766,25.2711034 c 164.34227,25.0749805 163.922398,25.0466617 163.542798,24.8773955 c 163.194478,24.7220773 163.083972,24.3552062 162.507536,24.2281877 c 161.9311,24.1011693 161.486505,24.2278936 161.023442,24.0729418 c 160.683926,23.9593313 160.613722,23.5600169 160.247214,23.4378478 c 159.788032,23.284787 159.142879,23.3079033 158.470122,23.3479115 c 157.944612,23.379163 157.40226,23.4207213 156.918545,23.3966683 c 156.540507,23.3778701 156.198284,23.318997 155.927908,23.1838091 c 155.706656,23.0731829 155.491061,22.6254308 154.912418,22.5487133 c 154.333776,22.4719957 153.815269,22.5408899 153.477743,22.3953734 c 153.162857,22.2596177 153.059763,21.8669612 152.90768,21.8289403 c 152.619885,21.7569915 152.344296,21.7089837 152.079711,21.6777671 c 151.723489,21.6357386 151.387213,21.6241467 151.067951,21.6255422 c 150.539712,21.6278512 150.05805,21.6657143 149.609675,21.6600961 c 149.164256,21.6545149 148.751688,21.606024 148.358942,21.4371413 c 147.975987,21.272469 147.898826,21.0950527 147.20865,20.4176183 c 146.996476,20.2093614 146.799175,20.0264466 146.620482,19.8623133 c 146.217909,19.4925413 145.909784,19.2180944 145.738826,18.9639568 c 145.491984,18.5970134 145.583581,18.1453902 145.368182,17.7925597 c 145.152783,17.4397292 145.165926,17.7984953 144.896094,16.9889995 c 144.810593,16.7324952 145.287599,16.6014434 145.368182,16.2988237 c 145.639372,15.2804067 145.643957,15.1850966 146.165982,14.6594246 c 146.360127,14.4639228 146.490629,13.4976765 146.72609,13.2622164 c 146.958739,13.0295674 147.222439,12.8230921 147.459269,12.5905543 c 147.691908,12.3621314 147.898621,12.1085601 148.024506,11.7803283 c 148.278545,11.1179515 147.982167,11.0791053 148.358942,10.9053089 c 148.735717,10.7315125 148.771309,11.160942 149.221661,11.3992716 c 149.504567,11.5489873 149.822314,11.5203687 150.197941,11.5545166 c 150.54641,11.5861956 151.129297,11.5778787 151.716397,11.5676281 c 152.606024,11.5520956 153.505322,11.5321234 153.61334,11.6401409 c 153.766508,11.7933085 153.971577,12.0864654 154.206096,12.2037249 c 154.698315,12.4498344 155.30673,12.4266608 155.920073,12.3855933 c 156.614755,12.3390797 157.315759,12.2696115 157.86142,12.5424421 c 157.984728,12.6040962 158.233373,12.9583493 158.341468,13.0123969 c 159.024284,13.3538048 159.278204,13.2615788 159.639685,13.4588758 c 159.942308,13.6240473 160.113559,13.8965364 160.486478,14.2710927 c 160.686801,14.4722946 160.886652,14.685981 161.1847,14.9840285 c 161.365783,15.1651113 161.695765,15.5370768 162.058957,15.9326287 c 162.26739,16.1596336 162.486761,16.3944067 162.695202,16.6053267 c 162.98153,16.8950608 163.247234,17.1397848 163.435633,17.2575338 c 164.000162,17.6103643 164.021282,17.325205 164.339387,17.4115026 c 164.810238,17.5392378 165.264047,17.4414623 165.657123,17.1984805 c 165.99567,16.9892063 166.289166,16.6722189 166.509697,16.2988237 c 166.701958,15.973295 166.74103,15.4880551 166.734104,14.9139927 c 166.728084,14.4149938 166.68731,13.8488822 166.682182,13.2622164 c 166.679036,12.9023848 166.695197,12.5006898 166.687695,12.1294325 c 166.677432,11.6215405 166.622882,11.170613 166.414031,10.9617617 c 166.128028,10.6757592 166.159076,10.7905709 165.962407,10.3972329 c 165.856212,10.1848444 165.924877,9.59242206 165.454331,9.16938297 c 164.983785,8.74634388 164.974481,8.36180797 164.452292,8.26613685 c 164.066307,8.19542002 163.685223,8.14909852 163.306397,8.12011209 c 162.791973,8.08075023 162.281711,8.07335435 161.768994,8.08024528 c 161.296026,8.08660196 160.820969,8.10511607 160.338626,8.12190977 c 159.851417,8.13887291 159.356775,8.15408076 158.849345,8.15323116 c 158.30471,8.15231927 157.5235,8.19792475 157.02874,8.07423488 c 156.478673,7.93671794 156.410757,7.53795352 155.927908,7.3770041 c 155.448624,7.21724255 154.659909,7.3770041 153.7078,7.27887313 c 153.241169,7.23077885 153.177424,6.91365194 152.90768,6.7284576 c 152.637688,6.54309314 152.23055,6.41794086 151.566924,6.47441966 c 151.12078,6.51238933 150.83207,6.57321218 150.371953,7.12362771 c 149.802771,7.36355241 149.279692,7.3710543 149.012431,7.50468461 c 148.879089,7.94916158 148.434746,8.09334816 147.953018,8.15070193 c 147.364023,8.22082663 146.719139,8.16114322 146.518475,8.36180796 c 146.394821,8.48546169 146.314956,8.65382226 146.165982,8.80279659 c 145.948678,9.02010087 145.440191,8.89079158 145.130308,9.04573276 c 144.88043,9.17067194 144.66026,9.44863746 144.456107,9.65279072 c 144.346843,9.76205412 143.798172,9.67084492 143.304335,10.0297316 c 142.810498,10.3886184 142.157865,11.0432512 142.08985,11.1355582 c 141.957628,11.3150011 141.643424,11.5108743 141.5439,11.7099224 c 141.102752,12.5922176 141.357822,12.4703809 141.236493,12.6980144 c 141.176405,12.8107505 141.061358,13.0475396 140.650022,13.4588759 c 140.642454,13.4664439 140.475326,13.8671445 140.479441,14.2710924 c 140.483232,14.6433395 140.496244,15.0055832 140.409424,15.1792235 c 140.115401,15.7672686 140.079292,15.5288997 139.848007,15.9935344 c 139.757941,16.1744701 139.58036,16.6889278 139.699324,16.9889992 c 139.885854,17.4594979 140.374476,17.9028886 140.409424,18.0776274 c 140.47668,18.413908 140.479334,18.662833 140.479441,19.1170398 c 140.479569,19.6665979 140.739019,19.9621383 141.05337,20.2266245 c 141.344686,20.4717288 141.683151,20.6901635 141.906073,21.0593792 c 142.369547,21.8270067 142.205884,22.3311986 142.369547,22.3953734 c 142.760603,22.5487133 143.328866,22.4966556 143.600649,22.768438 c 143.883821,23.05161 143.832385,23.9995403 144.093089,24.0729418 c 144.78673,24.2682369 144.91845,24.3326604 145.266426,24.4485295 c 145.614403,24.5643987 145.688025,24.9264879 146.165982,25.3320262 c 146.643938,25.7375646 146.888584,26.0064531 147.20865,26.3748412 c 147.473576,26.6797634 147.82896,26.7022626 148.178513,26.7926507 c 148.453283,26.863701 148.724449,26.9766989 148.945247,27.3017888 c 149.105281,27.5374137 149.59224,27.5897622 150.174708,27.5910411 c 150.625689,27.5920313 151.133927,27.5624062 151.592007,27.5635302 c 152.05168,27.5646582 152.460847,27.5967494 152.71097,27.7218108 c 153.119776,27.9262138 152.94786,27.9045245 153.190169,28.1468334 c 153.503682,28.460346 154.073659,28.4569839 154.654062,28.4357353 c 155.185456,28.416281 155.72559,28.3818333 156.085636,28.5618566 c 156.437966,28.7380211 156.369805,28.9284363 156.940166,29.2136168 c 157.074324,29.280696 157.497264,29.267697 157.908848,29.2925183 c 158.186808,29.3092812 158.459588,29.3432935 158.634741,29.4308696 c 159.069248,29.6481224 159.03255,29.8598897 159.195756,29.9142919 c 159.579797,30.0423054 160.386616,30.1669654 160.609554,30.3899032 c 160.99896,30.7793095 161.359188,31.2508628 161.647319,31.5165013 c 161.987582,31.8302031 162.21238,32.061759 162.332537,32.513503 c 162.403703,32.7810623 162.438162,33.1258629 162.438162,33.6107207 c 162.438162,34.0186619 162.452837,34.7374085 162.226661,34.9635845 c 162.099415,35.0908296 161.34304,35.1740848 161.184699,35.2532555 c 160.803153,35.4440287 160.436636,35.8990815 160.126547,36.2091704 c 159.841794,36.4939233 159.179177,37.2442251 158.849344,37.5271727 c 158.562323,37.7733947 157.94583,37.6887831 157.649859,37.8458129 c 157.288559,38.0375035 157.304334,38.2294232 157.028739,38.3672203 c 156.779615,38.4917827 156.349232,38.5325445 155.869915,38.5444681 c 155.482441,38.554107 155.062988,38.5449005 154.68146,38.5458838 c 154.13189,38.5473002 153.66101,38.5698593 153.477742,38.7003405 c 153.089172,38.9769904 153.019202,39.246292 152.448725,39.3350424 c 152.181539,39.3766091 151.804564,39.3785708 151.233648,39.3231334 c 150.418741,39.2440039 150.596106,38.8596605 150.197941,38.7003405 c 150.029005,38.632743 149.680384,38.5953741 149.321033,38.5608198 c 148.833433,38.5139334 148.326079,38.4722292 148.22107,38.3672203 c 147.962814,38.1089646 148.045898,38.0767299 147.786563,37.9548375 c 147.475527,37.8086441 146.977945,37.7125517 146.606729,37.4677074 c 146.218796,37.2118372 145.84714,36.8106482 145.465167,36.4029632 c 145.176871,36.0952613 144.882699,35.7838589 144.571216,35.5284429 c 144.330387,35.3309627 143.73117,35.1690929 143.304335,35.0360027 c 142.8775,34.9029125 143.0809,34.0577784 142.873936,33.8154271 c 142.600667,33.4954338 142.148359,33.3790233 141.446919,33.4138815 c 141.04774,33.4337188 140.567878,33.5025458 139.994415,33.6107207 c 139.62487,33.6804296 139.617843,34.0945714 138.951599,34.5884691 l 138.811517,35.1061011 z ";
// var path = point1 + point2 + point3 + point4 + point5 + point6 + point7 + point8 + point9;
var path = lettresANT + lettreS;

module.exports = path;
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/utilities.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/vector.js":[function(require,module,exports){
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
},{}],"/Users/Romain/Documents/Programmation/ants/AntColony/start.js":[function(require,module,exports){
'use strict';

var _antColony = require('./index.js');

var container = document.querySelector('.colony');

var options = {
	velocity: 0.001,
	nbAnts: 4000,
	weight: 10,
	repSize: 0.05,
	repSpeed: 0.002,
	nbStart: 500,
	nbRand: 900
	// obj par defaut
};

var antColony = _antColony(container, options);

window.addEventListener('click', function (){
	// options.velocity = 0.003;
	options.nbAnts = 20000;
	// options.weight = 10000000;
	// options.repSpeed = 0.01;
	// options.repSize = 0.1;

	// antColony.changeOptions(options);
	antColony.changeOptions(options);
});


},{"./index.js":"/Users/Romain/Documents/Programmation/ants/AntColony/index.js"}]},{},["/Users/Romain/Documents/Programmation/ants/AntColony/start.js"])
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi91c3IvbG9jYWwvbGliL25vZGVfbW9kdWxlcy93YXRjaGlmeS9ub2RlX21vZHVsZXMvYnJvd3NlcmlmeS9ub2RlX21vZHVsZXMvYnJvd3Nlci1wYWNrL19wcmVsdWRlLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9pbmRleC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9pY2guanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL25vZGVfbW9kdWxlcy90d28tc3VtL3R3by1zdW0uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL3JvYnVzdC1zY2FsZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL25vZGVfbW9kdWxlcy9yb2J1c3Qtc3VidHJhY3Qvcm9idXN0LWRpZmYuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1bS9yb2J1c3Qtc3VtLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vbm9kZV9tb2R1bGVzL3R3by1wcm9kdWN0L3R3by1wcm9kdWN0LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vb3JpZW50YXRpb24uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3NpbXBsaWNpYWwtY29tcGxleC9ub2RlX21vZHVsZXMvYml0LXR3aWRkbGUvdHdpZGRsZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvc2ltcGxpY2lhbC1jb21wbGV4L25vZGVfbW9kdWxlcy91bmlvbi1maW5kL2luZGV4LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvdG9wb2xvZ3kuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvdW5pcS91bmlxLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvdHJpYW5ndWxhdGUuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9wYXJzZS1zdmctcGF0aC9pbmRleC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2FudC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2FudHNHcm91cC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2NyZWF0ZUVkZ2VzLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvZWRnZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2luaXRpYWxpemVQb2ludHMuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L3NyYy9tb3VzZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3BvaW50LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvcmVuZGVyaW5nLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvc3ZnUGF0aC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3V0aWxpdGllcy5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3ZlY3Rvci5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3RhcnQuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2Q0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM3YkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNoQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNqREE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzNKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDM0pBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNoQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDN0xBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzVNQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN6REE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdFZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3pEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOUpBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdkRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDbk5BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNyRkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3hLQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ25CQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNUQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN0UEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNoQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM1QkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EiLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbiBlKHQsbixyKXtmdW5jdGlvbiBzKG8sdSl7aWYoIW5bb10pe2lmKCF0W29dKXt2YXIgYT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2lmKCF1JiZhKXJldHVybiBhKG8sITApO2lmKGkpcmV0dXJuIGkobywhMCk7dmFyIGY9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitvK1wiJ1wiKTt0aHJvdyBmLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsZn12YXIgbD1uW29dPXtleHBvcnRzOnt9fTt0W29dWzBdLmNhbGwobC5leHBvcnRzLGZ1bmN0aW9uKGUpe3ZhciBuPXRbb11bMV1bZV07cmV0dXJuIHMobj9uOmUpfSxsLGwuZXhwb3J0cyxlLHQsbixyKX1yZXR1cm4gbltvXS5leHBvcnRzfXZhciBpPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7Zm9yKHZhciBvPTA7bzxyLmxlbmd0aDtvKyspcyhyW29dKTtyZXR1cm4gc30pIiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgaW5pdFJlbmRlcmluZyA9IHJlcXVpcmUoJy4vc3JjL3JlbmRlcmluZy5qcycpO1xudmFyIGluaXRpYWxpemVQb2ludHMgPSByZXF1aXJlKCcuL3NyYy9pbml0aWFsaXplUG9pbnRzLmpzJyk7XG52YXIgY3JlYXRlRWRnZXMgPSByZXF1aXJlKCcuL3NyYy9jcmVhdGVFZGdlcy5qcycpO1xuLy8gdmFyIGluaXRBbnRzID0gcmVxdWlyZSgnLi9zcmMvaW5pdGlhbGl6ZUFudHMnKTtcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbiBpbml0KGNvbnRhaW5lckVsZW1lbnQsIG9wdGlvbnMpe1xuXG5cdHZhciByZW5kZXIsIHBvaW50c0luZm9zLCBlZGdlcywgcG9wdWxhdGlvbiwgcG9pbnRzTWFwO1xuXG5cblx0ZnVuY3Rpb24gX2luaXQoY29udGFpbmVyRWxlbWVudCwgb3B0aW9ucyl7XG5cdFx0cG9pbnRzSW5mb3MgPSBpbml0aWFsaXplUG9pbnRzKG9wdGlvbnMubmJTdGFydCwgb3B0aW9ucy5uYlJhbmQpO1xuXHRcdGVkZ2VzID0gY3JlYXRlRWRnZXMocG9pbnRzSW5mb3MucG9pbnRzKTtcblx0XHQvLyBwb3B1bGF0aW9uID0gb3B0aW9ucy5uYkFudHM7XG5cdFx0Ly8gcG9wdWxhdGlvbiA9IGluaXRBbnRzKGNvbnRhaW5lckVsZW1lbnQsIHBvaW50c0luZm9zLCBvcHRpb25zKTtcblx0XHRwb2ludHNNYXAgPSB7XG5cdFx0XHRwb2ludHNJbmZvczogcG9pbnRzSW5mb3MsXG5cdFx0XHRlZGdlczogZWRnZXNcblx0XHRcdC8vIHBvcHVsYXRpb246IHBvcHVsYXRpb25cblx0XHR9O1xuXHRcdHJlbmRlciA9IGluaXRSZW5kZXJpbmcoY29udGFpbmVyRWxlbWVudCwgcG9pbnRzTWFwLCBvcHRpb25zKTtcblx0fVxuXG5cdF9pbml0KGNvbnRhaW5lckVsZW1lbnQsIG9wdGlvbnMpO1xuXG5cdHJldHVybiB7XG5cdFx0dG9nZ2xlUGxheVBhdXNlOiBmdW5jdGlvbigpeyByZW5kZXIudG9nZ2xlUGxheVBhdXNlKCkgfSxcblx0XHRjaGFuZ2VPcHRpb25zOiBmdW5jdGlvbihvcHRzKXtcblx0XHRcdHJlbmRlci5tb2RpZnlBbnRzKG9wdHMpO1xuXHRcdH0sXG5cdFx0cmVzZXQ6IGZ1bmN0aW9uKG9wdHMpe1xuXHRcdFx0cmVuZGVyLnJlc2V0KCk7XG5cblx0XHRcdFx0Ly8gcmVzZXQgZWxlbWVudHNcblx0XHRcdF9pbml0KGNvbnRhaW5lckVsZW1lbnQsIG9wdHMpO1xuXHRcdH1cblx0fTtcbn07IiwiXCJ1c2Ugc3RyaWN0XCJcblxuLy9IaWdoIGxldmVsIGlkZWE6XG4vLyAxLiBVc2UgQ2xhcmtzb24ncyBpbmNyZW1lbnRhbCBjb25zdHJ1Y3Rpb24gdG8gZmluZCBjb252ZXggaHVsbFxuLy8gMi4gUG9pbnQgbG9jYXRpb24gaW4gdHJpYW5ndWxhdGlvbiBieSBqdW1wIGFuZCB3YWxrXG5cbm1vZHVsZS5leHBvcnRzID0gaW5jcmVtZW50YWxDb252ZXhIdWxsXG5cbnZhciBvcmllbnQgPSByZXF1aXJlKFwicm9idXN0LW9yaWVudGF0aW9uXCIpXG52YXIgY29tcGFyZUNlbGwgPSByZXF1aXJlKFwic2ltcGxpY2lhbC1jb21wbGV4XCIpLmNvbXBhcmVDZWxsc1xuXG5mdW5jdGlvbiBjb21wYXJlSW50KGEsIGIpIHtcbiAgcmV0dXJuIGEgLSBiXG59XG5cbmZ1bmN0aW9uIFNpbXBsZXgodmVydGljZXMsIGFkamFjZW50LCBib3VuZGFyeSkge1xuICB0aGlzLnZlcnRpY2VzID0gdmVydGljZXNcbiAgdGhpcy5hZGphY2VudCA9IGFkamFjZW50XG4gIHRoaXMuYm91bmRhcnkgPSBib3VuZGFyeVxuICB0aGlzLmxhc3RWaXNpdGVkID0gLTFcbn1cblxuU2ltcGxleC5wcm90b3R5cGUuZmxpcCA9IGZ1bmN0aW9uKCkge1xuICB2YXIgdCA9IHRoaXMudmVydGljZXNbMF1cbiAgdGhpcy52ZXJ0aWNlc1swXSA9IHRoaXMudmVydGljZXNbMV1cbiAgdGhpcy52ZXJ0aWNlc1sxXSA9IHRcbiAgdmFyIHUgPSB0aGlzLmFkamFjZW50WzBdXG4gIHRoaXMuYWRqYWNlbnRbMF0gPSB0aGlzLmFkamFjZW50WzFdXG4gIHRoaXMuYWRqYWNlbnRbMV0gPSB1XG59XG5cbmZ1bmN0aW9uIEdsdWVGYWNldCh2ZXJ0aWNlcywgY2VsbCwgaW5kZXgpIHtcbiAgdGhpcy52ZXJ0aWNlcyA9IHZlcnRpY2VzXG4gIHRoaXMuY2VsbCA9IGNlbGxcbiAgdGhpcy5pbmRleCA9IGluZGV4XG59XG5cbmZ1bmN0aW9uIGNvbXBhcmVHbHVlKGEsIGIpIHtcbiAgcmV0dXJuIGNvbXBhcmVDZWxsKGEudmVydGljZXMsIGIudmVydGljZXMpXG59XG5cbmZ1bmN0aW9uIGJha2VPcmllbnQoZCkge1xuICB2YXIgY29kZSA9IFtcImZ1bmN0aW9uIG9yaWVudCgpe3ZhciB0dXBsZT10aGlzLnR1cGxlO3JldHVybiB0ZXN0KFwiXVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgaWYoaSA+IDApIHtcbiAgICAgIGNvZGUucHVzaChcIixcIilcbiAgICB9XG4gICAgY29kZS5wdXNoKFwidHVwbGVbXCIsIGksIFwiXVwiKVxuICB9XG4gIGNvZGUucHVzaChcIil9cmV0dXJuIG9yaWVudFwiKVxuICB2YXIgcHJvYyA9IG5ldyBGdW5jdGlvbihcInRlc3RcIiwgY29kZS5qb2luKFwiXCIpKVxuICB2YXIgdGVzdCA9IG9yaWVudFtkKzFdXG4gIGlmKCF0ZXN0KSB7XG4gICAgdGVzdCA9IG9yaWVudFxuICB9XG4gIHJldHVybiBwcm9jKHRlc3QpXG59XG5cbnZhciBCQUtFRCA9IFtdXG5cbmZ1bmN0aW9uIFRyaWFuZ3VsYXRpb24oZGltZW5zaW9uLCB2ZXJ0aWNlcywgc2ltcGxpY2VzKSB7XG4gIHRoaXMuZGltZW5zaW9uID0gZGltZW5zaW9uXG4gIHRoaXMudmVydGljZXMgPSB2ZXJ0aWNlc1xuICB0aGlzLnNpbXBsaWNlcyA9IHNpbXBsaWNlc1xuICB0aGlzLmludGVyaW9yID0gc2ltcGxpY2VzLmZpbHRlcihmdW5jdGlvbihjKSB7XG4gICAgcmV0dXJuICFjLmJvdW5kYXJ5XG4gIH0pXG5cbiAgdGhpcy50dXBsZSA9IG5ldyBBcnJheShkaW1lbnNpb24rMSlcbiAgZm9yKHZhciBpPTA7IGk8PWRpbWVuc2lvbjsgKytpKSB7XG4gICAgdGhpcy50dXBsZVtpXSA9IHRoaXMudmVydGljZXNbaV1cbiAgfVxuXG4gIHZhciBvID0gQkFLRURbZGltZW5zaW9uXVxuICBpZighbykge1xuICAgIG8gPSBCQUtFRFtkaW1lbnNpb25dID0gYmFrZU9yaWVudChkaW1lbnNpb24pXG4gIH1cbiAgdGhpcy5vcmllbnQgPSBvXG59XG5cbnZhciBwcm90byA9IFRyaWFuZ3VsYXRpb24ucHJvdG90eXBlXG5cbi8vRGVnZW5lcmF0ZSBzaXR1YXRpb24gd2hlcmUgd2UgYXJlIG9uIGJvdW5kYXJ5LCBidXQgY29wbGFuYXIgdG8gZmFjZVxucHJvdG8uaGFuZGxlQm91bmRhcnlEZWdlbmVyYWN5ID0gZnVuY3Rpb24oY2VsbCwgcG9pbnQpIHtcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgbiA9IHRoaXMudmVydGljZXMubGVuZ3RoIC0gMVxuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcblxuICAvL0R1bWIgc29sdXRpb246IEp1c3QgZG8gZGZzIGZyb20gYm91bmRhcnkgY2VsbCB1bnRpbCB3ZSBmaW5kIGFueSBwZWFrLCBvciB0ZXJtaW5hdGVcbiAgdmFyIHRvVmlzaXQgPSBbIGNlbGwgXVxuICBjZWxsLmxhc3RWaXNpdGVkID0gLW5cbiAgd2hpbGUodG9WaXNpdC5sZW5ndGggPiAwKSB7XG4gICAgY2VsbCA9IHRvVmlzaXQucG9wKClcbiAgICB2YXIgY2VsbFZlcnRzID0gY2VsbC52ZXJ0aWNlc1xuICAgIHZhciBjZWxsQWRqID0gY2VsbC5hZGphY2VudFxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIHZhciBuZWlnaGJvciA9IGNlbGxBZGpbaV1cbiAgICAgIGlmKCFuZWlnaGJvci5ib3VuZGFyeSB8fCBuZWlnaGJvci5sYXN0VmlzaXRlZCA8PSAtbikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgdmFyIG52ID0gbmVpZ2hib3IudmVydGljZXNcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgdmFyIHZ2ID0gbnZbal1cbiAgICAgICAgaWYodnYgPCAwKSB7XG4gICAgICAgICAgdHVwbGVbal0gPSBwb2ludFxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIHR1cGxlW2pdID0gdmVydHNbdnZdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHZhciBvID0gdGhpcy5vcmllbnQoKVxuICAgICAgaWYobyA+IDApIHtcbiAgICAgICAgcmV0dXJuIG5laWdoYm9yXG4gICAgICB9XG4gICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IC1uXG4gICAgICBpZihvID09PSAwKSB7XG4gICAgICAgIHRvVmlzaXQucHVzaChuZWlnaGJvcilcbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgcmV0dXJuIG51bGxcbn1cblxucHJvdG8ud2FsayA9IGZ1bmN0aW9uKHBvaW50LCByYW5kb20pIHtcbiAgLy9BbGlhcyBsb2NhbCBwcm9wZXJ0aWVzXG4gIHZhciBuID0gdGhpcy52ZXJ0aWNlcy5sZW5ndGggLSAxXG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG5cbiAgLy9Db21wdXRlIGluaXRpYWwganVtcCBjZWxsXG4gIHZhciBpbml0SW5kZXggPSByYW5kb20gPyAodGhpcy5pbnRlcmlvci5sZW5ndGggKiBNYXRoLnJhbmRvbSgpKXwwIDogKHRoaXMuaW50ZXJpb3IubGVuZ3RoLTEpXG4gIHZhciBjZWxsID0gdGhpcy5pbnRlcmlvclsgaW5pdEluZGV4IF1cblxuICAvL1N0YXJ0IHdhbGtpbmdcbm91dGVyTG9vcDpcbiAgd2hpbGUoIWNlbGwuYm91bmRhcnkpIHtcbiAgICB2YXIgY2VsbFZlcnRzID0gY2VsbC52ZXJ0aWNlc1xuICAgIHZhciBjZWxsQWRqID0gY2VsbC5hZGphY2VudFxuXG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgdHVwbGVbaV0gPSB2ZXJ0c1tjZWxsVmVydHNbaV1dXG4gICAgfVxuICAgIGNlbGwubGFzdFZpc2l0ZWQgPSBuXG5cbiAgICAvL0ZpbmQgZmFydGhlc3QgYWRqYWNlbnQgY2VsbFxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIHZhciBuZWlnaGJvciA9IGNlbGxBZGpbaV1cbiAgICAgIGlmKG5laWdoYm9yLmxhc3RWaXNpdGVkID49IG4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHZhciBwcmV2ID0gdHVwbGVbaV1cbiAgICAgIHR1cGxlW2ldID0gcG9pbnRcbiAgICAgIHZhciBvID0gdGhpcy5vcmllbnQoKVxuICAgICAgdHVwbGVbaV0gPSBwcmV2XG4gICAgICBpZihvIDwgMCkge1xuICAgICAgICBjZWxsID0gbmVpZ2hib3JcbiAgICAgICAgY29udGludWUgb3V0ZXJMb29wXG4gICAgICB9IGVsc2Uge1xuICAgICAgICBpZighbmVpZ2hib3IuYm91bmRhcnkpIHtcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IG5cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IC1uXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gICAgcmV0dXJuXG4gIH1cblxuICByZXR1cm4gY2VsbFxufVxuXG5wcm90by5hZGRQZWFrcyA9IGZ1bmN0aW9uKHBvaW50LCBjZWxsKSB7XG4gIHZhciBuID0gdGhpcy52ZXJ0aWNlcy5sZW5ndGggLSAxXG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG4gIHZhciBpbnRlcmlvciA9IHRoaXMuaW50ZXJpb3JcbiAgdmFyIHNpbXBsaWNlcyA9IHRoaXMuc2ltcGxpY2VzXG5cbiAgLy9XYWxraW5nIGZpbmlzaGVkIGF0IGJvdW5kYXJ5LCB0aW1lIHRvIGFkZCBwZWFrc1xuICB2YXIgdG92aXNpdCA9IFsgY2VsbCBdXG5cbiAgLy9TdHJldGNoIGluaXRpYWwgYm91bmRhcnkgY2VsbCBpbnRvIGEgcGVha1xuICBjZWxsLmxhc3RWaXNpdGVkID0gblxuICBjZWxsLnZlcnRpY2VzW2NlbGwudmVydGljZXMuaW5kZXhPZigtMSldID0gblxuICBjZWxsLmJvdW5kYXJ5ID0gZmFsc2VcbiAgaW50ZXJpb3IucHVzaChjZWxsKVxuXG4gIC8vUmVjb3JkIGEgbGlzdCBvZiBhbGwgbmV3IGJvdW5kYXJpZXMgY3JlYXRlZCBieSBhZGRlZCBwZWFrcyBzbyB3ZSBjYW4gZ2x1ZSB0aGVtIHRvZ2V0aGVyIHdoZW4gd2UgYXJlIGFsbCBkb25lXG4gIHZhciBnbHVlRmFjZXRzID0gW11cblxuICAvL0RvIGEgdHJhdmVyc2FsIG9mIHRoZSBib3VuZGFyeSB3YWxraW5nIG91dHdhcmQgZnJvbSBzdGFydGluZyBwZWFrXG4gIHdoaWxlKHRvdmlzaXQubGVuZ3RoID4gMCkge1xuICAgIC8vUG9wIG9mZiBwZWFrIGFuZCB3YWxrIG92ZXIgYWRqYWNlbnQgY2VsbHNcbiAgICB2YXIgY2VsbCA9IHRvdmlzaXQucG9wKClcbiAgICB2YXIgY2VsbFZlcnRzID0gY2VsbC52ZXJ0aWNlc1xuICAgIHZhciBjZWxsQWRqID0gY2VsbC5hZGphY2VudFxuICAgIHZhciBpbmRleE9mTiA9IGNlbGxWZXJ0cy5pbmRleE9mKG4pXG4gICAgaWYoaW5kZXhPZk4gPCAwKSB7XG4gICAgICBjb250aW51ZVxuICAgIH1cblxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIGlmKGkgPT09IGluZGV4T2ZOKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG5cbiAgICAgIC8vRm9yIGVhY2ggYm91bmRhcnkgbmVpZ2hib3Igb2YgdGhlIGNlbGxcbiAgICAgIHZhciBuZWlnaGJvciA9IGNlbGxBZGpbaV1cbiAgICAgIGlmKCFuZWlnaGJvci5ib3VuZGFyeSB8fCBuZWlnaGJvci5sYXN0VmlzaXRlZCA+PSBuKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG5cbiAgICAgIHZhciBudiA9IG5laWdoYm9yLnZlcnRpY2VzXG5cbiAgICAgIC8vVGVzdCBpZiBuZWlnaGJvciBpcyBhIHBlYWtcbiAgICAgIGlmKG5laWdoYm9yLmxhc3RWaXNpdGVkICE9PSAtbikgeyAgICAgIFxuICAgICAgICAvL0NvbXB1dGUgb3JpZW50YXRpb24gb2YgcCByZWxhdGl2ZSB0byBlYWNoIGJvdW5kYXJ5IHBlYWtcbiAgICAgICAgdmFyIGluZGV4T2ZOZWcxID0gMFxuICAgICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgICAgaWYobnZbal0gPCAwKSB7XG4gICAgICAgICAgICBpbmRleE9mTmVnMSA9IGpcbiAgICAgICAgICAgIHR1cGxlW2pdID0gcG9pbnRcbiAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdHVwbGVbal0gPSB2ZXJ0c1tudltqXV1cbiAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgdmFyIG8gPSB0aGlzLm9yaWVudCgpXG5cbiAgICAgICAgLy9UZXN0IGlmIG5laWdoYm9yIGNlbGwgaXMgYWxzbyBhIHBlYWtcbiAgICAgICAgaWYobyA+IDApIHtcbiAgICAgICAgICBudltpbmRleE9mTmVnMV0gPSBuXG4gICAgICAgICAgbmVpZ2hib3IuYm91bmRhcnkgPSBmYWxzZVxuICAgICAgICAgIGludGVyaW9yLnB1c2gobmVpZ2hib3IpXG4gICAgICAgICAgdG92aXNpdC5wdXNoKG5laWdoYm9yKVxuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gblxuICAgICAgICAgIGNvbnRpbnVlXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSAtblxuICAgICAgICB9XG4gICAgICB9XG5cbiAgICAgIHZhciBuYSA9IG5laWdoYm9yLmFkamFjZW50XG5cbiAgICAgIC8vT3RoZXJ3aXNlLCByZXBsYWNlIG5laWdoYm9yIHdpdGggbmV3IGZhY2VcbiAgICAgIHZhciB2dmVydHMgPSBjZWxsVmVydHMuc2xpY2UoKVxuICAgICAgdmFyIHZhZGogPSBjZWxsQWRqLnNsaWNlKClcbiAgICAgIHZhciBuY2VsbCA9IG5ldyBTaW1wbGV4KHZ2ZXJ0cywgdmFkaiwgdHJ1ZSlcbiAgICAgIHNpbXBsaWNlcy5wdXNoKG5jZWxsKVxuXG4gICAgICAvL0Nvbm5lY3QgdG8gbmVpZ2hib3JcbiAgICAgIHZhciBvcHBvc2l0ZSA9IG5hLmluZGV4T2YoY2VsbClcbiAgICAgIGlmKG9wcG9zaXRlIDwgMCkge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbmFbb3Bwb3NpdGVdID0gbmNlbGxcbiAgICAgIHZhZGpbaW5kZXhPZk5dID0gbmVpZ2hib3JcblxuICAgICAgLy9Db25uZWN0IHRvIGNlbGxcbiAgICAgIHZ2ZXJ0c1tpXSA9IC0xXG4gICAgICB2YWRqW2ldID0gY2VsbFxuICAgICAgY2VsbEFkaltpXSA9IG5jZWxsXG5cbiAgICAgIC8vRmxpcCBmYWNldFxuICAgICAgbmNlbGwuZmxpcCgpXG5cbiAgICAgIC8vQWRkIHRvIGdsdWUgbGlzdFxuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICB2YXIgdXUgPSB2dmVydHNbal1cbiAgICAgICAgaWYodXUgPCAwIHx8IHV1ID09PSBuKSB7XG4gICAgICAgICAgY29udGludWVcbiAgICAgICAgfVxuICAgICAgICB2YXIgbmZhY2UgPSBuZXcgQXJyYXkoZC0xKVxuICAgICAgICB2YXIgbnB0ciA9IDBcbiAgICAgICAgZm9yKHZhciBrPTA7IGs8PWQ7ICsraykge1xuICAgICAgICAgIHZhciB2diA9IHZ2ZXJ0c1trXVxuICAgICAgICAgIGlmKHZ2IDwgMCB8fCBrID09PSBqKSB7XG4gICAgICAgICAgICBjb250aW51ZVxuICAgICAgICAgIH1cbiAgICAgICAgICBuZmFjZVtucHRyKytdID0gdnZcbiAgICAgICAgfVxuICAgICAgICBnbHVlRmFjZXRzLnB1c2gobmV3IEdsdWVGYWNldChuZmFjZSwgbmNlbGwsIGopKVxuICAgICAgfVxuICAgIH1cbiAgfVxuXG4gIC8vR2x1ZSBib3VuZGFyeSBmYWNldHMgdG9nZXRoZXJcbiAgZ2x1ZUZhY2V0cy5zb3J0KGNvbXBhcmVHbHVlKVxuXG4gIGZvcih2YXIgaT0wOyBpKzE8Z2x1ZUZhY2V0cy5sZW5ndGg7IGkrPTIpIHtcbiAgICB2YXIgYSA9IGdsdWVGYWNldHNbaV1cbiAgICB2YXIgYiA9IGdsdWVGYWNldHNbaSsxXVxuICAgIHZhciBhaSA9IGEuaW5kZXhcbiAgICB2YXIgYmkgPSBiLmluZGV4XG4gICAgaWYoYWkgPCAwIHx8IGJpIDwgMCkge1xuICAgICAgY29udGludWVcbiAgICB9XG4gICAgYS5jZWxsLmFkamFjZW50W2EuaW5kZXhdID0gYi5jZWxsXG4gICAgYi5jZWxsLmFkamFjZW50W2IuaW5kZXhdID0gYS5jZWxsXG4gIH1cbn1cblxucHJvdG8uaW5zZXJ0ID0gZnVuY3Rpb24ocG9pbnQsIHJhbmRvbSkge1xuICAvL0FkZCBwb2ludFxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG4gIHZlcnRzLnB1c2gocG9pbnQpXG5cbiAgdmFyIGNlbGwgPSB0aGlzLndhbGsocG9pbnQsIHJhbmRvbSlcbiAgaWYoIWNlbGwpIHtcbiAgICByZXR1cm5cbiAgfVxuXG4gIC8vQWxpYXMgbG9jYWwgcHJvcGVydGllc1xuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcblxuICAvL0RlZ2VuZXJhdGUgY2FzZTogSWYgcG9pbnQgaXMgY29wbGFuYXIgdG8gY2VsbCwgdGhlbiB3YWxrIHVudGlsIHdlIGZpbmQgYSBub24tZGVnZW5lcmF0ZSBib3VuZGFyeVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHZ2ID0gY2VsbC52ZXJ0aWNlc1tpXVxuICAgIGlmKHZ2IDwgMCkge1xuICAgICAgdHVwbGVbaV0gPSBwb2ludFxuICAgIH0gZWxzZSB7XG4gICAgICB0dXBsZVtpXSA9IHZlcnRzW3Z2XVxuICAgIH1cbiAgfVxuICB2YXIgbyA9IHRoaXMub3JpZW50KHR1cGxlKVxuICBpZihvIDwgMCkge1xuICAgIHJldHVyblxuICB9IGVsc2UgaWYobyA9PT0gMCkge1xuICAgIGNlbGwgPSB0aGlzLmhhbmRsZUJvdW5kYXJ5RGVnZW5lcmFjeShjZWxsLCBwb2ludClcbiAgICBpZighY2VsbCkge1xuICAgICAgcmV0dXJuXG4gICAgfVxuICB9XG5cbiAgLy9BZGQgcGVha3NcbiAgdGhpcy5hZGRQZWFrcyhwb2ludCwgY2VsbClcbn1cblxuLy9FeHRyYWN0IGFsbCBib3VuZGFyeSBjZWxsc1xucHJvdG8uYm91bmRhcnkgPSBmdW5jdGlvbigpIHtcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgYm91bmRhcnkgPSBbXVxuICB2YXIgY2VsbHMgPSB0aGlzLnNpbXBsaWNlc1xuICB2YXIgbmMgPSBjZWxscy5sZW5ndGhcbiAgZm9yKHZhciBpPTA7IGk8bmM7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBpZihjLmJvdW5kYXJ5KSB7XG4gICAgICB2YXIgYmNlbGwgPSBuZXcgQXJyYXkoZClcbiAgICAgIHZhciBjdiA9IGMudmVydGljZXNcbiAgICAgIHZhciBwdHIgPSAwXG4gICAgICB2YXIgcGFyaXR5ID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICBpZihjdltqXSA+PSAwKSB7XG4gICAgICAgICAgYmNlbGxbcHRyKytdID0gY3Zbal1cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICBwYXJpdHkgPSBqJjFcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgaWYocGFyaXR5ID09PSAoZCYxKSkge1xuICAgICAgICB2YXIgdCA9IGJjZWxsWzBdXG4gICAgICAgIGJjZWxsWzBdID0gYmNlbGxbMV1cbiAgICAgICAgYmNlbGxbMV0gPSB0XG4gICAgICB9XG4gICAgICBib3VuZGFyeS5wdXNoKGJjZWxsKVxuICAgIH1cbiAgfVxuICByZXR1cm4gYm91bmRhcnlcbn1cblxuZnVuY3Rpb24gaW5jcmVtZW50YWxDb252ZXhIdWxsKHBvaW50cywgcmFuZG9tU2VhcmNoKSB7XG4gIHZhciBuID0gcG9pbnRzLmxlbmd0aFxuICBpZihuID09PSAwKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKFwiTXVzdCBoYXZlIGF0IGxlYXN0IGQrMSBwb2ludHNcIilcbiAgfVxuICB2YXIgZCA9IHBvaW50c1swXS5sZW5ndGhcbiAgaWYobiA8PSBkKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKFwiTXVzdCBpbnB1dCBhdCBsZWFzdCBkKzEgcG9pbnRzXCIpXG4gIH1cblxuICAvL0ZJWE1FOiBUaGlzIGNvdWxkIGJlIGRlZ2VuZXJhdGUsIGJ1dCBuZWVkIHRvIHNlbGVjdCBkKzEgbm9uLWNvcGxhbmFyIHBvaW50cyB0byBib290c3RyYXAgcHJvY2Vzc1xuICB2YXIgaW5pdGlhbFNpbXBsZXggPSBwb2ludHMuc2xpY2UoMCwgZCsxKVxuXG4gIC8vTWFrZSBzdXJlIGluaXRpYWwgc2ltcGxleCBpcyBwb3NpdGl2ZWx5IG9yaWVudGVkXG4gIHZhciBvID0gb3JpZW50LmFwcGx5KHZvaWQgMCwgaW5pdGlhbFNpbXBsZXgpXG4gIGlmKG8gPT09IDApIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoXCJJbnB1dCBub3QgaW4gZ2VuZXJhbCBwb3NpdGlvblwiKVxuICB9XG4gIHZhciBpbml0aWFsQ29vcmRzID0gbmV3IEFycmF5KGQrMSlcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIGluaXRpYWxDb29yZHNbaV0gPSBpXG4gIH1cbiAgaWYobyA8IDApIHtcbiAgICBpbml0aWFsQ29vcmRzWzBdID0gMVxuICAgIGluaXRpYWxDb29yZHNbMV0gPSAwXG4gIH1cblxuICAvL0NyZWF0ZSBpbml0aWFsIHRvcG9sb2dpY2FsIGluZGV4LCBnbHVlIHBvaW50ZXJzIHRvZ2V0aGVyIChraW5kIG9mIG1lc3N5KVxuICB2YXIgaW5pdGlhbENlbGwgPSBuZXcgU2ltcGxleChpbml0aWFsQ29vcmRzLCBuZXcgQXJyYXkoZCsxKSwgZmFsc2UpXG4gIHZhciBib3VuZGFyeSA9IGluaXRpYWxDZWxsLmFkamFjZW50XG4gIHZhciBsaXN0ID0gbmV3IEFycmF5KGQrMilcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB2ZXJ0cyA9IGluaXRpYWxDb29yZHMuc2xpY2UoKVxuICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgIGlmKGogPT09IGkpIHtcbiAgICAgICAgdmVydHNbal0gPSAtMVxuICAgICAgfVxuICAgIH1cbiAgICB2YXIgdCA9IHZlcnRzWzBdXG4gICAgdmVydHNbMF0gPSB2ZXJ0c1sxXVxuICAgIHZlcnRzWzFdID0gdFxuICAgIHZhciBjZWxsID0gbmV3IFNpbXBsZXgodmVydHMsIG5ldyBBcnJheShkKzEpLCB0cnVlKVxuICAgIGJvdW5kYXJ5W2ldID0gY2VsbFxuICAgIGxpc3RbaV0gPSBjZWxsXG4gIH1cbiAgbGlzdFtkKzFdID0gaW5pdGlhbENlbGxcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB2ZXJ0cyA9IGJvdW5kYXJ5W2ldLnZlcnRpY2VzXG4gICAgdmFyIGFkaiA9IGJvdW5kYXJ5W2ldLmFkamFjZW50XG4gICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgdmFyIHYgPSB2ZXJ0c1tqXVxuICAgICAgaWYodiA8IDApIHtcbiAgICAgICAgYWRqW2pdID0gaW5pdGlhbENlbGxcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGZvcih2YXIgaz0wOyBrPD1kOyArK2spIHtcbiAgICAgICAgaWYoYm91bmRhcnlba10udmVydGljZXMuaW5kZXhPZih2KSA8IDApIHtcbiAgICAgICAgICBhZGpbal0gPSBib3VuZGFyeVtrXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICB9XG5cbiAgLy9Jbml0aWFsaXplIHRyaWFuZ2xlc1xuICB2YXIgdHJpYW5nbGVzID0gbmV3IFRyaWFuZ3VsYXRpb24oZCwgaW5pdGlhbFNpbXBsZXgsIGxpc3QpXG5cbiAgLy9JbnNlcnQgcmVtYWluaW5nIHBvaW50c1xuICB2YXIgdXNlUmFuZG9tID0gISFyYW5kb21TZWFyY2hcbiAgZm9yKHZhciBpPWQrMTsgaTxuOyArK2kpIHtcbiAgICB0cmlhbmdsZXMuaW5zZXJ0KHBvaW50c1tpXSwgdXNlUmFuZG9tKVxuICB9XG4gIFxuICAvL0V4dHJhY3QgYm91bmRhcnkgY2VsbHNcbiAgcmV0dXJuIHRyaWFuZ2xlcy5ib3VuZGFyeSgpXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSBmYXN0VHdvU3VtXG5cbmZ1bmN0aW9uIGZhc3RUd29TdW0oYSwgYiwgcmVzdWx0KSB7XG5cdHZhciB4ID0gYSArIGJcblx0dmFyIGJ2ID0geCAtIGFcblx0dmFyIGF2ID0geCAtIGJ2XG5cdHZhciBiciA9IGIgLSBidlxuXHR2YXIgYXIgPSBhIC0gYXZcblx0aWYocmVzdWx0KSB7XG5cdFx0cmVzdWx0WzBdID0gYXIgKyBiclxuXHRcdHJlc3VsdFsxXSA9IHhcblx0XHRyZXR1cm4gcmVzdWx0XG5cdH1cblx0cmV0dXJuIFthciticiwgeF1cbn0iLCJcInVzZSBzdHJpY3RcIlxuXG52YXIgdHdvUHJvZHVjdCA9IHJlcXVpcmUoXCJ0d28tcHJvZHVjdFwiKVxudmFyIHR3b1N1bSA9IHJlcXVpcmUoXCJ0d28tc3VtXCIpXG5cbm1vZHVsZS5leHBvcnRzID0gc2NhbGVMaW5lYXJFeHBhbnNpb25cblxuZnVuY3Rpb24gc2NhbGVMaW5lYXJFeHBhbnNpb24oZSwgc2NhbGUpIHtcbiAgdmFyIG4gPSBlLmxlbmd0aFxuICBpZihuID09PSAxKSB7XG4gICAgdmFyIHRzID0gdHdvUHJvZHVjdChlWzBdLCBzY2FsZSlcbiAgICBpZih0c1swXSkge1xuICAgICAgcmV0dXJuIHRzXG4gICAgfVxuICAgIHJldHVybiBbIHRzWzFdIF1cbiAgfVxuICB2YXIgZyA9IG5ldyBBcnJheSgyICogbilcbiAgdmFyIHEgPSBbMC4xLCAwLjFdXG4gIHZhciB0ID0gWzAuMSwgMC4xXVxuICB2YXIgY291bnQgPSAwXG4gIHR3b1Byb2R1Y3QoZVswXSwgc2NhbGUsIHEpXG4gIGlmKHFbMF0pIHtcbiAgICBnW2NvdW50KytdID0gcVswXVxuICB9XG4gIGZvcih2YXIgaT0xOyBpPG47ICsraSkge1xuICAgIHR3b1Byb2R1Y3QoZVtpXSwgc2NhbGUsIHQpXG4gICAgdmFyIHBxID0gcVsxXVxuICAgIHR3b1N1bShwcSwgdFswXSwgcSlcbiAgICBpZihxWzBdKSB7XG4gICAgICBnW2NvdW50KytdID0gcVswXVxuICAgIH1cbiAgICB2YXIgYSA9IHRbMV1cbiAgICB2YXIgYiA9IHFbMV1cbiAgICB2YXIgeCA9IGEgKyBiXG4gICAgdmFyIGJ2ID0geCAtIGFcbiAgICB2YXIgeSA9IGIgLSBidlxuICAgIHFbMV0gPSB4XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gIH1cbiAgaWYocVsxXSkge1xuICAgIGdbY291bnQrK10gPSBxWzFdXG4gIH1cbiAgaWYoY291bnQgPT09IDApIHtcbiAgICBnW2NvdW50KytdID0gMC4wXG4gIH1cbiAgZy5sZW5ndGggPSBjb3VudFxuICByZXR1cm4gZ1xufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gcm9idXN0U3VidHJhY3RcblxuLy9FYXN5IGNhc2U6IEFkZCB0d28gc2NhbGFyc1xuZnVuY3Rpb24gc2NhbGFyU2NhbGFyKGEsIGIpIHtcbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgYXYgPSB4IC0gYnZcbiAgdmFyIGJyID0gYiAtIGJ2XG4gIHZhciBhciA9IGEgLSBhdlxuICB2YXIgeSA9IGFyICsgYnJcbiAgaWYoeSkge1xuICAgIHJldHVybiBbeSwgeF1cbiAgfVxuICByZXR1cm4gW3hdXG59XG5cbmZ1bmN0aW9uIHJvYnVzdFN1YnRyYWN0KGUsIGYpIHtcbiAgdmFyIG5lID0gZS5sZW5ndGh8MFxuICB2YXIgbmYgPSBmLmxlbmd0aHwwXG4gIGlmKG5lID09PSAxICYmIG5mID09PSAxKSB7XG4gICAgcmV0dXJuIHNjYWxhclNjYWxhcihlWzBdLCAtZlswXSlcbiAgfVxuICB2YXIgbiA9IG5lICsgbmZcbiAgdmFyIGcgPSBuZXcgQXJyYXkobilcbiAgdmFyIGNvdW50ID0gMFxuICB2YXIgZXB0ciA9IDBcbiAgdmFyIGZwdHIgPSAwXG4gIHZhciBhYnMgPSBNYXRoLmFic1xuICB2YXIgZWkgPSBlW2VwdHJdXG4gIHZhciBlYSA9IGFicyhlaSlcbiAgdmFyIGZpID0gLWZbZnB0cl1cbiAgdmFyIGZhID0gYWJzKGZpKVxuICB2YXIgYSwgYlxuICBpZihlYSA8IGZhKSB7XG4gICAgYiA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBiID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gLWZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICBpZigoZXB0ciA8IG5lICYmIGVhIDwgZmEpIHx8IChmcHRyID49IG5mKSkge1xuICAgIGEgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYSA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgeSA9IGIgLSBidlxuICB2YXIgcTAgPSB5XG4gIHZhciBxMSA9IHhcbiAgdmFyIF94LCBfYnYsIF9hdiwgX2JyLCBfYXJcbiAgd2hpbGUoZXB0ciA8IG5lICYmIGZwdHIgPCBuZikge1xuICAgIGlmKGVhIDwgZmEpIHtcbiAgICAgIGEgPSBlaVxuICAgICAgZXB0ciArPSAxXG4gICAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICAgIGVhID0gYWJzKGVpKVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICBhID0gZmlcbiAgICAgIGZwdHIgKz0gMVxuICAgICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICAgIGZpID0gLWZbZnB0cl1cbiAgICAgICAgZmEgPSBhYnMoZmkpXG4gICAgICB9XG4gICAgfVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgfVxuICB3aGlsZShlcHRyIDwgbmUpIHtcbiAgICBhID0gZWlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICB9XG4gIH1cbiAgd2hpbGUoZnB0ciA8IG5mKSB7XG4gICAgYSA9IGZpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH0gXG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gLWZbZnB0cl1cbiAgICB9XG4gIH1cbiAgaWYocTApIHtcbiAgICBnW2NvdW50KytdID0gcTBcbiAgfVxuICBpZihxMSkge1xuICAgIGdbY291bnQrK10gPSBxMVxuICB9XG4gIGlmKCFjb3VudCkge1xuICAgIGdbY291bnQrK10gPSAwLjAgIFxuICB9XG4gIGcubGVuZ3RoID0gY291bnRcbiAgcmV0dXJuIGdcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IGxpbmVhckV4cGFuc2lvblN1bVxuXG4vL0Vhc3kgY2FzZTogQWRkIHR3byBzY2FsYXJzXG5mdW5jdGlvbiBzY2FsYXJTY2FsYXIoYSwgYikge1xuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciBhdiA9IHggLSBidlxuICB2YXIgYnIgPSBiIC0gYnZcbiAgdmFyIGFyID0gYSAtIGF2XG4gIHZhciB5ID0gYXIgKyBiclxuICBpZih5KSB7XG4gICAgcmV0dXJuIFt5LCB4XVxuICB9XG4gIHJldHVybiBbeF1cbn1cblxuZnVuY3Rpb24gbGluZWFyRXhwYW5zaW9uU3VtKGUsIGYpIHtcbiAgdmFyIG5lID0gZS5sZW5ndGh8MFxuICB2YXIgbmYgPSBmLmxlbmd0aHwwXG4gIGlmKG5lID09PSAxICYmIG5mID09PSAxKSB7XG4gICAgcmV0dXJuIHNjYWxhclNjYWxhcihlWzBdLCBmWzBdKVxuICB9XG4gIHZhciBuID0gbmUgKyBuZlxuICB2YXIgZyA9IG5ldyBBcnJheShuKVxuICB2YXIgY291bnQgPSAwXG4gIHZhciBlcHRyID0gMFxuICB2YXIgZnB0ciA9IDBcbiAgdmFyIGFicyA9IE1hdGguYWJzXG4gIHZhciBlaSA9IGVbZXB0cl1cbiAgdmFyIGVhID0gYWJzKGVpKVxuICB2YXIgZmkgPSBmW2ZwdHJdXG4gIHZhciBmYSA9IGFicyhmaSlcbiAgdmFyIGEsIGJcbiAgaWYoZWEgPCBmYSkge1xuICAgIGIgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYiA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IGZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICBpZigoZXB0ciA8IG5lICYmIGVhIDwgZmEpIHx8IChmcHRyID49IG5mKSkge1xuICAgIGEgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYSA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IGZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciB5ID0gYiAtIGJ2XG4gIHZhciBxMCA9IHlcbiAgdmFyIHExID0geFxuICB2YXIgX3gsIF9idiwgX2F2LCBfYnIsIF9hclxuICB3aGlsZShlcHRyIDwgbmUgJiYgZnB0ciA8IG5mKSB7XG4gICAgaWYoZWEgPCBmYSkge1xuICAgICAgYSA9IGVpXG4gICAgICBlcHRyICs9IDFcbiAgICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgICBlaSA9IGVbZXB0cl1cbiAgICAgICAgZWEgPSBhYnMoZWkpXG4gICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgIGEgPSBmaVxuICAgICAgZnB0ciArPSAxXG4gICAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgICAgIGZhID0gYWJzKGZpKVxuICAgICAgfVxuICAgIH1cbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gIH1cbiAgd2hpbGUoZXB0ciA8IG5lKSB7XG4gICAgYSA9IGVpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgfVxuICB9XG4gIHdoaWxlKGZwdHIgPCBuZikge1xuICAgIGEgPSBmaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9IFxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IGZbZnB0cl1cbiAgICB9XG4gIH1cbiAgaWYocTApIHtcbiAgICBnW2NvdW50KytdID0gcTBcbiAgfVxuICBpZihxMSkge1xuICAgIGdbY291bnQrK10gPSBxMVxuICB9XG4gIGlmKCFjb3VudCkge1xuICAgIGdbY291bnQrK10gPSAwLjAgIFxuICB9XG4gIGcubGVuZ3RoID0gY291bnRcbiAgcmV0dXJuIGdcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IHR3b1Byb2R1Y3RcblxudmFyIFNQTElUVEVSID0gKyhNYXRoLnBvdygyLCAyNykgKyAxLjApXG5cbmZ1bmN0aW9uIHR3b1Byb2R1Y3QoYSwgYiwgcmVzdWx0KSB7XG4gIHZhciB4ID0gYSAqIGJcblxuICB2YXIgYyA9IFNQTElUVEVSICogYVxuICB2YXIgYWJpZyA9IGMgLSBhXG4gIHZhciBhaGkgPSBjIC0gYWJpZ1xuICB2YXIgYWxvID0gYSAtIGFoaVxuXG4gIHZhciBkID0gU1BMSVRURVIgKiBiXG4gIHZhciBiYmlnID0gZCAtIGJcbiAgdmFyIGJoaSA9IGQgLSBiYmlnXG4gIHZhciBibG8gPSBiIC0gYmhpXG5cbiAgdmFyIGVycjEgPSB4IC0gKGFoaSAqIGJoaSlcbiAgdmFyIGVycjIgPSBlcnIxIC0gKGFsbyAqIGJoaSlcbiAgdmFyIGVycjMgPSBlcnIyIC0gKGFoaSAqIGJsbylcblxuICB2YXIgeSA9IGFsbyAqIGJsbyAtIGVycjNcblxuICBpZihyZXN1bHQpIHtcbiAgICByZXN1bHRbMF0gPSB5XG4gICAgcmVzdWx0WzFdID0geFxuICAgIHJldHVybiByZXN1bHRcbiAgfVxuXG4gIHJldHVybiBbIHksIHggXVxufSIsIlwidXNlIHN0cmljdFwiXG5cbnZhciB0d29Qcm9kdWN0ID0gcmVxdWlyZShcInR3by1wcm9kdWN0XCIpXG52YXIgcm9idXN0U3VtID0gcmVxdWlyZShcInJvYnVzdC1zdW1cIilcbnZhciByb2J1c3RTY2FsZSA9IHJlcXVpcmUoXCJyb2J1c3Qtc2NhbGVcIilcbnZhciByb2J1c3RTdWJ0cmFjdCA9IHJlcXVpcmUoXCJyb2J1c3Qtc3VidHJhY3RcIilcblxudmFyIE5VTV9FWFBBTkQgPSA1XG5cbnZhciBFUFNJTE9OICAgICA9IDEuMTEwMjIzMDI0NjI1MTU2NWUtMTZcbnZhciBFUlJCT1VORDMgICA9ICgzLjAgKyAxNi4wICogRVBTSUxPTikgKiBFUFNJTE9OXG52YXIgRVJSQk9VTkQ0ICAgPSAoNy4wICsgNTYuMCAqIEVQU0lMT04pICogRVBTSUxPTlxuXG5mdW5jdGlvbiBjb2ZhY3RvcihtLCBjKSB7XG4gIHZhciByZXN1bHQgPSBuZXcgQXJyYXkobS5sZW5ndGgtMSlcbiAgZm9yKHZhciBpPTE7IGk8bS5sZW5ndGg7ICsraSkge1xuICAgIHZhciByID0gcmVzdWx0W2ktMV0gPSBuZXcgQXJyYXkobS5sZW5ndGgtMSlcbiAgICBmb3IodmFyIGo9MCxrPTA7IGo8bS5sZW5ndGg7ICsraikge1xuICAgICAgaWYoaiA9PT0gYykge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgcltrKytdID0gbVtpXVtqXVxuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzdWx0XG59XG5cbmZ1bmN0aW9uIG1hdHJpeChuKSB7XG4gIHZhciByZXN1bHQgPSBuZXcgQXJyYXkobilcbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgcmVzdWx0W2ldID0gbmV3IEFycmF5KG4pXG4gICAgZm9yKHZhciBqPTA7IGo8bjsgKytqKSB7XG4gICAgICByZXN1bHRbaV1bal0gPSBbXCJtXCIsIGosIFwiW1wiLCAobi1pLTEpLCBcIl1cIl0uam9pbihcIlwiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzdWx0XG59XG5cbmZ1bmN0aW9uIHNpZ24obikge1xuICBpZihuICYgMSkge1xuICAgIHJldHVybiBcIi1cIlxuICB9XG4gIHJldHVybiBcIlwiXG59XG5cbmZ1bmN0aW9uIGdlbmVyYXRlU3VtKGV4cHIpIHtcbiAgaWYoZXhwci5sZW5ndGggPT09IDEpIHtcbiAgICByZXR1cm4gZXhwclswXVxuICB9IGVsc2UgaWYoZXhwci5sZW5ndGggPT09IDIpIHtcbiAgICByZXR1cm4gW1wic3VtKFwiLCBleHByWzBdLCBcIixcIiwgZXhwclsxXSwgXCIpXCJdLmpvaW4oXCJcIilcbiAgfSBlbHNlIHtcbiAgICB2YXIgbSA9IGV4cHIubGVuZ3RoPj4xXG4gICAgcmV0dXJuIFtcInN1bShcIiwgZ2VuZXJhdGVTdW0oZXhwci5zbGljZSgwLCBtKSksIFwiLFwiLCBnZW5lcmF0ZVN1bShleHByLnNsaWNlKG0pKSwgXCIpXCJdLmpvaW4oXCJcIilcbiAgfVxufVxuXG5mdW5jdGlvbiBkZXRlcm1pbmFudChtKSB7XG4gIGlmKG0ubGVuZ3RoID09PSAyKSB7XG4gICAgcmV0dXJuIFtbXCJzdW0ocHJvZChcIiwgbVswXVswXSwgXCIsXCIsIG1bMV1bMV0sIFwiKSxwcm9kKC1cIiwgbVswXVsxXSwgXCIsXCIsIG1bMV1bMF0sIFwiKSlcIl0uam9pbihcIlwiKV1cbiAgfSBlbHNlIHtcbiAgICB2YXIgZXhwciA9IFtdXG4gICAgZm9yKHZhciBpPTA7IGk8bS5sZW5ndGg7ICsraSkge1xuICAgICAgZXhwci5wdXNoKFtcInNjYWxlKFwiLCBnZW5lcmF0ZVN1bShkZXRlcm1pbmFudChjb2ZhY3RvcihtLCBpKSkpLCBcIixcIiwgc2lnbihpKSwgbVswXVtpXSwgXCIpXCJdLmpvaW4oXCJcIikpXG4gICAgfVxuICAgIHJldHVybiBleHByXG4gIH1cbn1cblxuZnVuY3Rpb24gb3JpZW50YXRpb24obikge1xuICB2YXIgcG9zID0gW11cbiAgdmFyIG5lZyA9IFtdXG4gIHZhciBtID0gbWF0cml4KG4pXG4gIHZhciBhcmdzID0gW11cbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgaWYoKGkmMSk9PT0wKSB7XG4gICAgICBwb3MucHVzaC5hcHBseShwb3MsIGRldGVybWluYW50KGNvZmFjdG9yKG0sIGkpKSlcbiAgICB9IGVsc2Uge1xuICAgICAgbmVnLnB1c2guYXBwbHkobmVnLCBkZXRlcm1pbmFudChjb2ZhY3RvcihtLCBpKSkpXG4gICAgfVxuICAgIGFyZ3MucHVzaChcIm1cIiArIGkpXG4gIH1cbiAgdmFyIHBvc0V4cHIgPSBnZW5lcmF0ZVN1bShwb3MpXG4gIHZhciBuZWdFeHByID0gZ2VuZXJhdGVTdW0obmVnKVxuICB2YXIgZnVuY05hbWUgPSBcIm9yaWVudGF0aW9uXCIgKyBuICsgXCJFeGFjdFwiXG4gIHZhciBjb2RlID0gW1wiZnVuY3Rpb24gXCIsIGZ1bmNOYW1lLCBcIihcIiwgYXJncy5qb2luKCksIFwiKXt2YXIgcD1cIiwgcG9zRXhwciwgXCIsbj1cIiwgbmVnRXhwciwgXCIsZD1zdWIocCxuKTtcXFxucmV0dXJuIGRbZC5sZW5ndGgtMV07fTtyZXR1cm4gXCIsIGZ1bmNOYW1lXS5qb2luKFwiXCIpXG4gIHZhciBwcm9jID0gbmV3IEZ1bmN0aW9uKFwic3VtXCIsIFwicHJvZFwiLCBcInNjYWxlXCIsIFwic3ViXCIsIGNvZGUpXG4gIHJldHVybiBwcm9jKHJvYnVzdFN1bSwgdHdvUHJvZHVjdCwgcm9idXN0U2NhbGUsIHJvYnVzdFN1YnRyYWN0KVxufVxuXG52YXIgb3JpZW50YXRpb24zRXhhY3QgPSBvcmllbnRhdGlvbigzKVxudmFyIG9yaWVudGF0aW9uNEV4YWN0ID0gb3JpZW50YXRpb24oNClcblxudmFyIENBQ0hFRCA9IFtcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24wKCkgeyByZXR1cm4gMCB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjEoKSB7IHJldHVybiAwIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMihhLCBiKSB7IFxuICAgIHJldHVybiBiWzBdIC0gYVswXVxuICB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjMoYSwgYiwgYykge1xuICAgIHZhciBsID0gKGFbMV0gLSBjWzFdKSAqIChiWzBdIC0gY1swXSlcbiAgICB2YXIgciA9IChhWzBdIC0gY1swXSkgKiAoYlsxXSAtIGNbMV0pXG4gICAgdmFyIGRldCA9IGwgLSByXG4gICAgdmFyIHNcbiAgICBpZihsID4gMCkge1xuICAgICAgaWYociA8PSAwKSB7XG4gICAgICAgIHJldHVybiBkZXRcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIHMgPSBsICsgclxuICAgICAgfVxuICAgIH0gZWxzZSBpZihsIDwgMCkge1xuICAgICAgaWYociA+PSAwKSB7XG4gICAgICAgIHJldHVybiBkZXRcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIHMgPSAtKGwgKyByKVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICByZXR1cm4gZGV0XG4gICAgfVxuICAgIHZhciB0b2wgPSBFUlJCT1VORDMgKiBzXG4gICAgaWYoZGV0ID49IHRvbCB8fCBkZXQgPD0gLXRvbCkge1xuICAgICAgcmV0dXJuIGRldFxuICAgIH1cbiAgICByZXR1cm4gb3JpZW50YXRpb24zRXhhY3QoYSwgYiwgYylcbiAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb240KGEsYixjLGQpIHtcbiAgICB2YXIgYWR4ID0gYVswXSAtIGRbMF1cbiAgICB2YXIgYmR4ID0gYlswXSAtIGRbMF1cbiAgICB2YXIgY2R4ID0gY1swXSAtIGRbMF1cbiAgICB2YXIgYWR5ID0gYVsxXSAtIGRbMV1cbiAgICB2YXIgYmR5ID0gYlsxXSAtIGRbMV1cbiAgICB2YXIgY2R5ID0gY1sxXSAtIGRbMV1cbiAgICB2YXIgYWR6ID0gYVsyXSAtIGRbMl1cbiAgICB2YXIgYmR6ID0gYlsyXSAtIGRbMl1cbiAgICB2YXIgY2R6ID0gY1syXSAtIGRbMl1cbiAgICB2YXIgYmR4Y2R5ID0gYmR4ICogY2R5XG4gICAgdmFyIGNkeGJkeSA9IGNkeCAqIGJkeVxuICAgIHZhciBjZHhhZHkgPSBjZHggKiBhZHlcbiAgICB2YXIgYWR4Y2R5ID0gYWR4ICogY2R5XG4gICAgdmFyIGFkeGJkeSA9IGFkeCAqIGJkeVxuICAgIHZhciBiZHhhZHkgPSBiZHggKiBhZHlcbiAgICB2YXIgZGV0ID0gYWR6ICogKGJkeGNkeSAtIGNkeGJkeSkgXG4gICAgICAgICAgICArIGJkeiAqIChjZHhhZHkgLSBhZHhjZHkpXG4gICAgICAgICAgICArIGNkeiAqIChhZHhiZHkgLSBiZHhhZHkpXG4gICAgdmFyIHBlcm1hbmVudCA9IChNYXRoLmFicyhiZHhjZHkpICsgTWF0aC5hYnMoY2R4YmR5KSkgKiBNYXRoLmFicyhhZHopXG4gICAgICAgICAgICAgICAgICArIChNYXRoLmFicyhjZHhhZHkpICsgTWF0aC5hYnMoYWR4Y2R5KSkgKiBNYXRoLmFicyhiZHopXG4gICAgICAgICAgICAgICAgICArIChNYXRoLmFicyhhZHhiZHkpICsgTWF0aC5hYnMoYmR4YWR5KSkgKiBNYXRoLmFicyhjZHopXG4gICAgdmFyIHRvbCA9IEVSUkJPVU5ENCAqIHBlcm1hbmVudFxuICAgIGlmICgoZGV0ID4gdG9sKSB8fCAoLWRldCA+IHRvbCkpIHtcbiAgICAgIHJldHVybiBkZXRcbiAgICB9XG4gICAgcmV0dXJuIG9yaWVudGF0aW9uNEV4YWN0KGEsYixjLGQpXG4gIH1cbl1cblxuZnVuY3Rpb24gc2xvd09yaWVudChhcmdzKSB7XG4gIHZhciBwcm9jID0gQ0FDSEVEW2FyZ3MubGVuZ3RoXVxuICBpZighcHJvYykge1xuICAgIHByb2MgPSBDQUNIRURbYXJncy5sZW5ndGhdID0gb3JpZW50YXRpb24oYXJncy5sZW5ndGgpXG4gIH1cbiAgcmV0dXJuIHByb2MuYXBwbHkodW5kZWZpbmVkLCBhcmdzKVxufVxuXG5mdW5jdGlvbiBnZW5lcmF0ZU9yaWVudGF0aW9uUHJvYygpIHtcbiAgd2hpbGUoQ0FDSEVELmxlbmd0aCA8PSBOVU1fRVhQQU5EKSB7XG4gICAgQ0FDSEVELnB1c2gob3JpZW50YXRpb24oQ0FDSEVELmxlbmd0aCkpXG4gIH1cbiAgdmFyIGFyZ3MgPSBbXVxuICB2YXIgcHJvY0FyZ3MgPSBbXCJzbG93XCJdXG4gIGZvcih2YXIgaT0wOyBpPD1OVU1fRVhQQU5EOyArK2kpIHtcbiAgICBhcmdzLnB1c2goXCJhXCIgKyBpKVxuICAgIHByb2NBcmdzLnB1c2goXCJvXCIgKyBpKVxuICB9XG4gIHZhciBjb2RlID0gW1xuICAgIFwiZnVuY3Rpb24gZ2V0T3JpZW50YXRpb24oXCIsIGFyZ3Muam9pbigpLCBcIil7c3dpdGNoKGFyZ3VtZW50cy5sZW5ndGgpe2Nhc2UgMDpjYXNlIDE6cmV0dXJuIDA7XCJcbiAgXVxuICBmb3IodmFyIGk9MjsgaTw9TlVNX0VYUEFORDsgKytpKSB7XG4gICAgY29kZS5wdXNoKFwiY2FzZSBcIiwgaSwgXCI6cmV0dXJuIG9cIiwgaSwgXCIoXCIsIGFyZ3Muc2xpY2UoMCwgaSkuam9pbigpLCBcIik7XCIpXG4gIH1cbiAgY29kZS5wdXNoKFwifXZhciBzPW5ldyBBcnJheShhcmd1bWVudHMubGVuZ3RoKTtmb3IodmFyIGk9MDtpPGFyZ3VtZW50cy5sZW5ndGg7KytpKXtzW2ldPWFyZ3VtZW50c1tpXX07cmV0dXJuIHNsb3cocyk7fXJldHVybiBnZXRPcmllbnRhdGlvblwiKVxuICBwcm9jQXJncy5wdXNoKGNvZGUuam9pbihcIlwiKSlcblxuICB2YXIgcHJvYyA9IEZ1bmN0aW9uLmFwcGx5KHVuZGVmaW5lZCwgcHJvY0FyZ3MpXG4gIG1vZHVsZS5leHBvcnRzID0gcHJvYy5hcHBseSh1bmRlZmluZWQsIFtzbG93T3JpZW50XS5jb25jYXQoQ0FDSEVEKSlcbiAgZm9yKHZhciBpPTA7IGk8PU5VTV9FWFBBTkQ7ICsraSkge1xuICAgIG1vZHVsZS5leHBvcnRzW2ldID0gQ0FDSEVEW2ldXG4gIH1cbn1cblxuZ2VuZXJhdGVPcmllbnRhdGlvblByb2MoKSIsIi8qKlxuICogQml0IHR3aWRkbGluZyBoYWNrcyBmb3IgSmF2YVNjcmlwdC5cbiAqXG4gKiBBdXRob3I6IE1pa29sYSBMeXNlbmtvXG4gKlxuICogUG9ydGVkIGZyb20gU3RhbmZvcmQgYml0IHR3aWRkbGluZyBoYWNrIGxpYnJhcnk6XG4gKiAgICBodHRwOi8vZ3JhcGhpY3Muc3RhbmZvcmQuZWR1L35zZWFuZGVyL2JpdGhhY2tzLmh0bWxcbiAqL1xuXG5cInVzZSBzdHJpY3RcIjsgXCJ1c2UgcmVzdHJpY3RcIjtcblxuLy9OdW1iZXIgb2YgYml0cyBpbiBhbiBpbnRlZ2VyXG52YXIgSU5UX0JJVFMgPSAzMjtcblxuLy9Db25zdGFudHNcbmV4cG9ydHMuSU5UX0JJVFMgID0gSU5UX0JJVFM7XG5leHBvcnRzLklOVF9NQVggICA9ICAweDdmZmZmZmZmO1xuZXhwb3J0cy5JTlRfTUlOICAgPSAtMTw8KElOVF9CSVRTLTEpO1xuXG4vL1JldHVybnMgLTEsIDAsICsxIGRlcGVuZGluZyBvbiBzaWduIG9mIHhcbmV4cG9ydHMuc2lnbiA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICh2ID4gMCkgLSAodiA8IDApO1xufVxuXG4vL0NvbXB1dGVzIGFic29sdXRlIHZhbHVlIG9mIGludGVnZXJcbmV4cG9ydHMuYWJzID0gZnVuY3Rpb24odikge1xuICB2YXIgbWFzayA9IHYgPj4gKElOVF9CSVRTLTEpO1xuICByZXR1cm4gKHYgXiBtYXNrKSAtIG1hc2s7XG59XG5cbi8vQ29tcHV0ZXMgbWluaW11bSBvZiBpbnRlZ2VycyB4IGFuZCB5XG5leHBvcnRzLm1pbiA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgcmV0dXJuIHkgXiAoKHggXiB5KSAmIC0oeCA8IHkpKTtcbn1cblxuLy9Db21wdXRlcyBtYXhpbXVtIG9mIGludGVnZXJzIHggYW5kIHlcbmV4cG9ydHMubWF4ID0gZnVuY3Rpb24oeCwgeSkge1xuICByZXR1cm4geCBeICgoeCBeIHkpICYgLSh4IDwgeSkpO1xufVxuXG4vL0NoZWNrcyBpZiBhIG51bWJlciBpcyBhIHBvd2VyIG9mIHR3b1xuZXhwb3J0cy5pc1BvdzIgPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAhKHYgJiAodi0xKSkgJiYgKCEhdik7XG59XG5cbi8vQ29tcHV0ZXMgbG9nIGJhc2UgMiBvZiB2XG5leHBvcnRzLmxvZzIgPSBmdW5jdGlvbih2KSB7XG4gIHZhciByLCBzaGlmdDtcbiAgciA9ICAgICAodiA+IDB4RkZGRikgPDwgNDsgdiA+Pj49IHI7XG4gIHNoaWZ0ID0gKHYgPiAweEZGICApIDw8IDM7IHYgPj4+PSBzaGlmdDsgciB8PSBzaGlmdDtcbiAgc2hpZnQgPSAodiA+IDB4RiAgICkgPDwgMjsgdiA+Pj49IHNoaWZ0OyByIHw9IHNoaWZ0O1xuICBzaGlmdCA9ICh2ID4gMHgzICAgKSA8PCAxOyB2ID4+Pj0gc2hpZnQ7IHIgfD0gc2hpZnQ7XG4gIHJldHVybiByIHwgKHYgPj4gMSk7XG59XG5cbi8vQ29tcHV0ZXMgbG9nIGJhc2UgMTAgb2YgdlxuZXhwb3J0cy5sb2cxMCA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICAodiA+PSAxMDAwMDAwMDAwKSA/IDkgOiAodiA+PSAxMDAwMDAwMDApID8gOCA6ICh2ID49IDEwMDAwMDAwKSA/IDcgOlxuICAgICAgICAgICh2ID49IDEwMDAwMDApID8gNiA6ICh2ID49IDEwMDAwMCkgPyA1IDogKHYgPj0gMTAwMDApID8gNCA6XG4gICAgICAgICAgKHYgPj0gMTAwMCkgPyAzIDogKHYgPj0gMTAwKSA/IDIgOiAodiA+PSAxMCkgPyAxIDogMDtcbn1cblxuLy9Db3VudHMgbnVtYmVyIG9mIGJpdHNcbmV4cG9ydHMucG9wQ291bnQgPSBmdW5jdGlvbih2KSB7XG4gIHYgPSB2IC0gKCh2ID4+PiAxKSAmIDB4NTU1NTU1NTUpO1xuICB2ID0gKHYgJiAweDMzMzMzMzMzKSArICgodiA+Pj4gMikgJiAweDMzMzMzMzMzKTtcbiAgcmV0dXJuICgodiArICh2ID4+PiA0KSAmIDB4RjBGMEYwRikgKiAweDEwMTAxMDEpID4+PiAyNDtcbn1cblxuLy9Db3VudHMgbnVtYmVyIG9mIHRyYWlsaW5nIHplcm9zXG5mdW5jdGlvbiBjb3VudFRyYWlsaW5nWmVyb3Modikge1xuICB2YXIgYyA9IDMyO1xuICB2ICY9IC12O1xuICBpZiAodikgYy0tO1xuICBpZiAodiAmIDB4MDAwMEZGRkYpIGMgLT0gMTY7XG4gIGlmICh2ICYgMHgwMEZGMDBGRikgYyAtPSA4O1xuICBpZiAodiAmIDB4MEYwRjBGMEYpIGMgLT0gNDtcbiAgaWYgKHYgJiAweDMzMzMzMzMzKSBjIC09IDI7XG4gIGlmICh2ICYgMHg1NTU1NTU1NSkgYyAtPSAxO1xuICByZXR1cm4gYztcbn1cbmV4cG9ydHMuY291bnRUcmFpbGluZ1plcm9zID0gY291bnRUcmFpbGluZ1plcm9zO1xuXG4vL1JvdW5kcyB0byBuZXh0IHBvd2VyIG9mIDJcbmV4cG9ydHMubmV4dFBvdzIgPSBmdW5jdGlvbih2KSB7XG4gIHYgKz0gdiA9PT0gMDtcbiAgLS12O1xuICB2IHw9IHYgPj4+IDE7XG4gIHYgfD0gdiA+Pj4gMjtcbiAgdiB8PSB2ID4+PiA0O1xuICB2IHw9IHYgPj4+IDg7XG4gIHYgfD0gdiA+Pj4gMTY7XG4gIHJldHVybiB2ICsgMTtcbn1cblxuLy9Sb3VuZHMgZG93biB0byBwcmV2aW91cyBwb3dlciBvZiAyXG5leHBvcnRzLnByZXZQb3cyID0gZnVuY3Rpb24odikge1xuICB2IHw9IHYgPj4+IDE7XG4gIHYgfD0gdiA+Pj4gMjtcbiAgdiB8PSB2ID4+PiA0O1xuICB2IHw9IHYgPj4+IDg7XG4gIHYgfD0gdiA+Pj4gMTY7XG4gIHJldHVybiB2IC0gKHY+Pj4xKTtcbn1cblxuLy9Db21wdXRlcyBwYXJpdHkgb2Ygd29yZFxuZXhwb3J0cy5wYXJpdHkgPSBmdW5jdGlvbih2KSB7XG4gIHYgXj0gdiA+Pj4gMTY7XG4gIHYgXj0gdiA+Pj4gODtcbiAgdiBePSB2ID4+PiA0O1xuICB2ICY9IDB4ZjtcbiAgcmV0dXJuICgweDY5OTYgPj4+IHYpICYgMTtcbn1cblxudmFyIFJFVkVSU0VfVEFCTEUgPSBuZXcgQXJyYXkoMjU2KTtcblxuKGZ1bmN0aW9uKHRhYikge1xuICBmb3IodmFyIGk9MDsgaTwyNTY7ICsraSkge1xuICAgIHZhciB2ID0gaSwgciA9IGksIHMgPSA3O1xuICAgIGZvciAodiA+Pj49IDE7IHY7IHYgPj4+PSAxKSB7XG4gICAgICByIDw8PSAxO1xuICAgICAgciB8PSB2ICYgMTtcbiAgICAgIC0tcztcbiAgICB9XG4gICAgdGFiW2ldID0gKHIgPDwgcykgJiAweGZmO1xuICB9XG59KShSRVZFUlNFX1RBQkxFKTtcblxuLy9SZXZlcnNlIGJpdHMgaW4gYSAzMiBiaXQgd29yZFxuZXhwb3J0cy5yZXZlcnNlID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gIChSRVZFUlNFX1RBQkxFWyB2ICAgICAgICAgJiAweGZmXSA8PCAyNCkgfFxuICAgICAgICAgIChSRVZFUlNFX1RBQkxFWyh2ID4+PiA4KSAgJiAweGZmXSA8PCAxNikgfFxuICAgICAgICAgIChSRVZFUlNFX1RBQkxFWyh2ID4+PiAxNikgJiAweGZmXSA8PCA4KSAgfFxuICAgICAgICAgICBSRVZFUlNFX1RBQkxFWyh2ID4+PiAyNCkgJiAweGZmXTtcbn1cblxuLy9JbnRlcmxlYXZlIGJpdHMgb2YgMiBjb29yZGluYXRlcyB3aXRoIDE2IGJpdHMuICBVc2VmdWwgZm9yIGZhc3QgcXVhZHRyZWUgY29kZXNcbmV4cG9ydHMuaW50ZXJsZWF2ZTIgPSBmdW5jdGlvbih4LCB5KSB7XG4gIHggJj0gMHhGRkZGO1xuICB4ID0gKHggfCAoeCA8PCA4KSkgJiAweDAwRkYwMEZGO1xuICB4ID0gKHggfCAoeCA8PCA0KSkgJiAweDBGMEYwRjBGO1xuICB4ID0gKHggfCAoeCA8PCAyKSkgJiAweDMzMzMzMzMzO1xuICB4ID0gKHggfCAoeCA8PCAxKSkgJiAweDU1NTU1NTU1O1xuXG4gIHkgJj0gMHhGRkZGO1xuICB5ID0gKHkgfCAoeSA8PCA4KSkgJiAweDAwRkYwMEZGO1xuICB5ID0gKHkgfCAoeSA8PCA0KSkgJiAweDBGMEYwRjBGO1xuICB5ID0gKHkgfCAoeSA8PCAyKSkgJiAweDMzMzMzMzMzO1xuICB5ID0gKHkgfCAoeSA8PCAxKSkgJiAweDU1NTU1NTU1O1xuXG4gIHJldHVybiB4IHwgKHkgPDwgMSk7XG59XG5cbi8vRXh0cmFjdHMgdGhlIG50aCBpbnRlcmxlYXZlZCBjb21wb25lbnRcbmV4cG9ydHMuZGVpbnRlcmxlYXZlMiA9IGZ1bmN0aW9uKHYsIG4pIHtcbiAgdiA9ICh2ID4+PiBuKSAmIDB4NTU1NTU1NTU7XG4gIHYgPSAodiB8ICh2ID4+PiAxKSkgICYgMHgzMzMzMzMzMztcbiAgdiA9ICh2IHwgKHYgPj4+IDIpKSAgJiAweDBGMEYwRjBGO1xuICB2ID0gKHYgfCAodiA+Pj4gNCkpICAmIDB4MDBGRjAwRkY7XG4gIHYgPSAodiB8ICh2ID4+PiAxNikpICYgMHgwMDBGRkZGO1xuICByZXR1cm4gKHYgPDwgMTYpID4+IDE2O1xufVxuXG5cbi8vSW50ZXJsZWF2ZSBiaXRzIG9mIDMgY29vcmRpbmF0ZXMsIGVhY2ggd2l0aCAxMCBiaXRzLiAgVXNlZnVsIGZvciBmYXN0IG9jdHJlZSBjb2Rlc1xuZXhwb3J0cy5pbnRlcmxlYXZlMyA9IGZ1bmN0aW9uKHgsIHksIHopIHtcbiAgeCAmPSAweDNGRjtcbiAgeCAgPSAoeCB8ICh4PDwxNikpICYgNDI3ODE5MDMzNTtcbiAgeCAgPSAoeCB8ICh4PDw4KSkgICYgMjUxNzE5Njk1O1xuICB4ICA9ICh4IHwgKHg8PDQpKSAgJiAzMjcyMzU2MDM1O1xuICB4ICA9ICh4IHwgKHg8PDIpKSAgJiAxMjI3MTMzNTEzO1xuXG4gIHkgJj0gMHgzRkY7XG4gIHkgID0gKHkgfCAoeTw8MTYpKSAmIDQyNzgxOTAzMzU7XG4gIHkgID0gKHkgfCAoeTw8OCkpICAmIDI1MTcxOTY5NTtcbiAgeSAgPSAoeSB8ICh5PDw0KSkgICYgMzI3MjM1NjAzNTtcbiAgeSAgPSAoeSB8ICh5PDwyKSkgICYgMTIyNzEzMzUxMztcbiAgeCB8PSAoeSA8PCAxKTtcbiAgXG4gIHogJj0gMHgzRkY7XG4gIHogID0gKHogfCAoejw8MTYpKSAmIDQyNzgxOTAzMzU7XG4gIHogID0gKHogfCAoejw8OCkpICAmIDI1MTcxOTY5NTtcbiAgeiAgPSAoeiB8ICh6PDw0KSkgICYgMzI3MjM1NjAzNTtcbiAgeiAgPSAoeiB8ICh6PDwyKSkgICYgMTIyNzEzMzUxMztcbiAgXG4gIHJldHVybiB4IHwgKHogPDwgMik7XG59XG5cbi8vRXh0cmFjdHMgbnRoIGludGVybGVhdmVkIGNvbXBvbmVudCBvZiBhIDMtdHVwbGVcbmV4cG9ydHMuZGVpbnRlcmxlYXZlMyA9IGZ1bmN0aW9uKHYsIG4pIHtcbiAgdiA9ICh2ID4+PiBuKSAgICAgICAmIDEyMjcxMzM1MTM7XG4gIHYgPSAodiB8ICh2Pj4+MikpICAgJiAzMjcyMzU2MDM1O1xuICB2ID0gKHYgfCAodj4+PjQpKSAgICYgMjUxNzE5Njk1O1xuICB2ID0gKHYgfCAodj4+PjgpKSAgICYgNDI3ODE5MDMzNTtcbiAgdiA9ICh2IHwgKHY+Pj4xNikpICAmIDB4M0ZGO1xuICByZXR1cm4gKHY8PDIyKT4+MjI7XG59XG5cbi8vQ29tcHV0ZXMgbmV4dCBjb21iaW5hdGlvbiBpbiBjb2xleGljb2dyYXBoaWMgb3JkZXIgKHRoaXMgaXMgbWlzdGFrZW5seSBjYWxsZWQgbmV4dFBlcm11dGF0aW9uIG9uIHRoZSBiaXQgdHdpZGRsaW5nIGhhY2tzIHBhZ2UpXG5leHBvcnRzLm5leHRDb21iaW5hdGlvbiA9IGZ1bmN0aW9uKHYpIHtcbiAgdmFyIHQgPSB2IHwgKHYgLSAxKTtcbiAgcmV0dXJuICh0ICsgMSkgfCAoKCh+dCAmIC1+dCkgLSAxKSA+Pj4gKGNvdW50VHJhaWxpbmdaZXJvcyh2KSArIDEpKTtcbn1cblxuIiwiXCJ1c2Ugc3RyaWN0XCI7IFwidXNlIHJlc3RyaWN0XCI7XG5cbm1vZHVsZS5leHBvcnRzID0gVW5pb25GaW5kO1xuXG5mdW5jdGlvbiBVbmlvbkZpbmQoY291bnQpIHtcbiAgdGhpcy5yb290cyA9IG5ldyBBcnJheShjb3VudCk7XG4gIHRoaXMucmFua3MgPSBuZXcgQXJyYXkoY291bnQpO1xuICBcbiAgZm9yKHZhciBpPTA7IGk8Y291bnQ7ICsraSkge1xuICAgIHRoaXMucm9vdHNbaV0gPSBpO1xuICAgIHRoaXMucmFua3NbaV0gPSAwO1xuICB9XG59XG5cbnZhciBwcm90byA9IFVuaW9uRmluZC5wcm90b3R5cGVcblxuT2JqZWN0LmRlZmluZVByb3BlcnR5KHByb3RvLCBcImxlbmd0aFwiLCB7XG4gIFwiZ2V0XCI6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnJvb3RzLmxlbmd0aFxuICB9XG59KVxuXG5wcm90by5tYWtlU2V0ID0gZnVuY3Rpb24oKSB7XG4gIHZhciBuID0gdGhpcy5yb290cy5sZW5ndGg7XG4gIHRoaXMucm9vdHMucHVzaChuKTtcbiAgdGhpcy5yYW5rcy5wdXNoKDApO1xuICByZXR1cm4gbjtcbn1cblxucHJvdG8uZmluZCA9IGZ1bmN0aW9uKHgpIHtcbiAgdmFyIHJvb3RzID0gdGhpcy5yb290cztcbiAgd2hpbGUocm9vdHNbeF0gIT09IHgpIHtcbiAgICB2YXIgeSA9IHJvb3RzW3hdO1xuICAgIHJvb3RzW3hdID0gcm9vdHNbeV07XG4gICAgeCA9IHk7XG4gIH1cbiAgcmV0dXJuIHg7XG59XG5cbnByb3RvLmxpbmsgPSBmdW5jdGlvbih4LCB5KSB7XG4gIHZhciB4ciA9IHRoaXMuZmluZCh4KVxuICAgICwgeXIgPSB0aGlzLmZpbmQoeSk7XG4gIGlmKHhyID09PSB5cikge1xuICAgIHJldHVybjtcbiAgfVxuICB2YXIgcmFua3MgPSB0aGlzLnJhbmtzXG4gICAgLCByb290cyA9IHRoaXMucm9vdHNcbiAgICAsIHhkICAgID0gcmFua3NbeHJdXG4gICAgLCB5ZCAgICA9IHJhbmtzW3lyXTtcbiAgaWYoeGQgPCB5ZCkge1xuICAgIHJvb3RzW3hyXSA9IHlyO1xuICB9IGVsc2UgaWYoeWQgPCB4ZCkge1xuICAgIHJvb3RzW3lyXSA9IHhyO1xuICB9IGVsc2Uge1xuICAgIHJvb3RzW3lyXSA9IHhyO1xuICAgICsrcmFua3NbeHJdO1xuICB9XG59IiwiXCJ1c2Ugc3RyaWN0XCI7IFwidXNlIHJlc3RyaWN0XCI7XG5cbnZhciBiaXRzICAgICAgPSByZXF1aXJlKFwiYml0LXR3aWRkbGVcIilcbiAgLCBVbmlvbkZpbmQgPSByZXF1aXJlKFwidW5pb24tZmluZFwiKVxuXG4vL1JldHVybnMgdGhlIGRpbWVuc2lvbiBvZiBhIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gZGltZW5zaW9uKGNlbGxzKSB7XG4gIHZhciBkID0gMFxuICAgICwgbWF4ID0gTWF0aC5tYXhcbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgZCA9IG1heChkLCBjZWxsc1tpXS5sZW5ndGgpXG4gIH1cbiAgcmV0dXJuIGQtMVxufVxuZXhwb3J0cy5kaW1lbnNpb24gPSBkaW1lbnNpb25cblxuLy9Db3VudHMgdGhlIG51bWJlciBvZiB2ZXJ0aWNlcyBpbiBmYWNlc1xuZnVuY3Rpb24gY291bnRWZXJ0aWNlcyhjZWxscykge1xuICB2YXIgdmMgPSAtMVxuICAgICwgbWF4ID0gTWF0aC5tYXhcbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wLCBqbD1jLmxlbmd0aDsgajxqbDsgKytqKSB7XG4gICAgICB2YyA9IG1heCh2YywgY1tqXSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHZjKzFcbn1cbmV4cG9ydHMuY291bnRWZXJ0aWNlcyA9IGNvdW50VmVydGljZXNcblxuLy9SZXR1cm5zIGEgZGVlcCBjb3B5IG9mIGNlbGxzXG5mdW5jdGlvbiBjbG9uZUNlbGxzKGNlbGxzKSB7XG4gIHZhciBuY2VsbHMgPSBuZXcgQXJyYXkoY2VsbHMubGVuZ3RoKVxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICBuY2VsbHNbaV0gPSBjZWxsc1tpXS5zbGljZSgwKVxuICB9XG4gIHJldHVybiBuY2VsbHNcbn1cbmV4cG9ydHMuY2xvbmVDZWxscyA9IGNsb25lQ2VsbHNcblxuLy9SYW5rcyBhIHBhaXIgb2YgY2VsbHMgdXAgdG8gcGVybXV0YXRpb25cbmZ1bmN0aW9uIGNvbXBhcmVDZWxscyhhLCBiKSB7XG4gIHZhciBuID0gYS5sZW5ndGhcbiAgICAsIHQgPSBhLmxlbmd0aCAtIGIubGVuZ3RoXG4gICAgLCBtaW4gPSBNYXRoLm1pblxuICBpZih0KSB7XG4gICAgcmV0dXJuIHRcbiAgfVxuICBzd2l0Y2gobikge1xuICAgIGNhc2UgMDpcbiAgICAgIHJldHVybiAwO1xuICAgIGNhc2UgMTpcbiAgICAgIHJldHVybiBhWzBdIC0gYlswXTtcbiAgICBjYXNlIDI6XG4gICAgICB2YXIgZCA9IGFbMF0rYVsxXS1iWzBdLWJbMV1cbiAgICAgIGlmKGQpIHtcbiAgICAgICAgcmV0dXJuIGRcbiAgICAgIH1cbiAgICAgIHJldHVybiBtaW4oYVswXSxhWzFdKSAtIG1pbihiWzBdLGJbMV0pXG4gICAgY2FzZSAzOlxuICAgICAgdmFyIGwxID0gYVswXSthWzFdXG4gICAgICAgICwgbTEgPSBiWzBdK2JbMV1cbiAgICAgIGQgPSBsMSthWzJdIC0gKG0xK2JbMl0pXG4gICAgICBpZihkKSB7XG4gICAgICAgIHJldHVybiBkXG4gICAgICB9XG4gICAgICB2YXIgbDAgPSBtaW4oYVswXSwgYVsxXSlcbiAgICAgICAgLCBtMCA9IG1pbihiWzBdLCBiWzFdKVxuICAgICAgICAsIGQgID0gbWluKGwwLCBhWzJdKSAtIG1pbihtMCwgYlsyXSlcbiAgICAgIGlmKGQpIHtcbiAgICAgICAgcmV0dXJuIGRcbiAgICAgIH1cbiAgICAgIHJldHVybiBtaW4obDArYVsyXSwgbDEpIC0gbWluKG0wK2JbMl0sIG0xKVxuICAgIFxuICAgIC8vVE9ETzogTWF5YmUgb3B0aW1pemUgbj00IGFzIHdlbGw/XG4gICAgXG4gICAgZGVmYXVsdDpcbiAgICAgIHZhciBhcyA9IGEuc2xpY2UoMClcbiAgICAgIGFzLnNvcnQoKVxuICAgICAgdmFyIGJzID0gYi5zbGljZSgwKVxuICAgICAgYnMuc29ydCgpXG4gICAgICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICAgICAgdCA9IGFzW2ldIC0gYnNbaV1cbiAgICAgICAgaWYodCkge1xuICAgICAgICAgIHJldHVybiB0XG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJldHVybiAwXG4gIH1cbn1cbmV4cG9ydHMuY29tcGFyZUNlbGxzID0gY29tcGFyZUNlbGxzXG5cbmZ1bmN0aW9uIGNvbXBhcmVaaXBwZWQoYSwgYikge1xuICByZXR1cm4gY29tcGFyZUNlbGxzKGFbMF0sIGJbMF0pXG59XG5cbi8vUHV0cyBhIGNlbGwgY29tcGxleCBpbnRvIG5vcm1hbCBvcmRlciBmb3IgdGhlIHB1cnBvc2VzIG9mIGZpbmRDZWxsIHF1ZXJpZXNcbmZ1bmN0aW9uIG5vcm1hbGl6ZShjZWxscywgYXR0cikge1xuICBpZihhdHRyKSB7XG4gICAgdmFyIGxlbiA9IGNlbGxzLmxlbmd0aFxuICAgIHZhciB6aXBwZWQgPSBuZXcgQXJyYXkobGVuKVxuICAgIGZvcih2YXIgaT0wOyBpPGxlbjsgKytpKSB7XG4gICAgICB6aXBwZWRbaV0gPSBbY2VsbHNbaV0sIGF0dHJbaV1dXG4gICAgfVxuICAgIHppcHBlZC5zb3J0KGNvbXBhcmVaaXBwZWQpXG4gICAgZm9yKHZhciBpPTA7IGk8bGVuOyArK2kpIHtcbiAgICAgIGNlbGxzW2ldID0gemlwcGVkW2ldWzBdXG4gICAgICBhdHRyW2ldID0gemlwcGVkW2ldWzFdXG4gICAgfVxuICAgIHJldHVybiBjZWxsc1xuICB9IGVsc2Uge1xuICAgIGNlbGxzLnNvcnQoY29tcGFyZUNlbGxzKVxuICAgIHJldHVybiBjZWxsc1xuICB9XG59XG5leHBvcnRzLm5vcm1hbGl6ZSA9IG5vcm1hbGl6ZVxuXG4vL1JlbW92ZXMgYWxsIGR1cGxpY2F0ZSBjZWxscyBpbiB0aGUgY29tcGxleFxuZnVuY3Rpb24gdW5pcXVlKGNlbGxzKSB7XG4gIGlmKGNlbGxzLmxlbmd0aCA9PT0gMCkge1xuICAgIHJldHVybiBbXVxuICB9XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBjZWxscy5sZW5ndGhcbiAgZm9yKHZhciBpPTE7IGk8bGVuOyArK2kpIHtcbiAgICB2YXIgYSA9IGNlbGxzW2ldXG4gICAgaWYoY29tcGFyZUNlbGxzKGEsIGNlbGxzW2ktMV0pKSB7XG4gICAgICBpZihpID09PSBwdHIpIHtcbiAgICAgICAgcHRyKytcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGNlbGxzW3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgY2VsbHMubGVuZ3RoID0gcHRyXG4gIHJldHVybiBjZWxsc1xufVxuZXhwb3J0cy51bmlxdWUgPSB1bmlxdWU7XG5cbi8vRmluZHMgYSBjZWxsIGluIGEgbm9ybWFsaXplZCBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGZpbmRDZWxsKGNlbGxzLCBjKSB7XG4gIHZhciBsbyA9IDBcbiAgICAsIGhpID0gY2VsbHMubGVuZ3RoLTFcbiAgICAsIHIgID0gLTFcbiAgd2hpbGUgKGxvIDw9IGhpKSB7XG4gICAgdmFyIG1pZCA9IChsbyArIGhpKSA+PiAxXG4gICAgICAsIHMgICA9IGNvbXBhcmVDZWxscyhjZWxsc1ttaWRdLCBjKVxuICAgIGlmKHMgPD0gMCkge1xuICAgICAgaWYocyA9PT0gMCkge1xuICAgICAgICByID0gbWlkXG4gICAgICB9XG4gICAgICBsbyA9IG1pZCArIDFcbiAgICB9IGVsc2UgaWYocyA+IDApIHtcbiAgICAgIGhpID0gbWlkIC0gMVxuICAgIH1cbiAgfVxuICByZXR1cm4gclxufVxuZXhwb3J0cy5maW5kQ2VsbCA9IGZpbmRDZWxsO1xuXG4vL0J1aWxkcyBhbiBpbmRleCBmb3IgYW4gbi1jZWxsLiAgVGhpcyBpcyBtb3JlIGdlbmVyYWwgdGhhbiBkdWFsLCBidXQgbGVzcyBlZmZpY2llbnRcbmZ1bmN0aW9uIGluY2lkZW5jZShmcm9tX2NlbGxzLCB0b19jZWxscykge1xuICB2YXIgaW5kZXggPSBuZXcgQXJyYXkoZnJvbV9jZWxscy5sZW5ndGgpXG4gIGZvcih2YXIgaT0wLCBpbD1pbmRleC5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIGluZGV4W2ldID0gW11cbiAgfVxuICB2YXIgYiA9IFtdXG4gIGZvcih2YXIgaT0wLCBuPXRvX2NlbGxzLmxlbmd0aDsgaTxuOyArK2kpIHtcbiAgICB2YXIgYyA9IHRvX2NlbGxzW2ldXG4gICAgdmFyIGNsID0gYy5sZW5ndGhcbiAgICBmb3IodmFyIGs9MSwga249KDE8PGNsKTsgazxrbjsgKytrKSB7XG4gICAgICBiLmxlbmd0aCA9IGJpdHMucG9wQ291bnQoaylcbiAgICAgIHZhciBsID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8Y2w7ICsraikge1xuICAgICAgICBpZihrICYgKDE8PGopKSB7XG4gICAgICAgICAgYltsKytdID0gY1tqXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICB2YXIgaWR4PWZpbmRDZWxsKGZyb21fY2VsbHMsIGIpXG4gICAgICBpZihpZHggPCAwKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICB3aGlsZSh0cnVlKSB7XG4gICAgICAgIGluZGV4W2lkeCsrXS5wdXNoKGkpXG4gICAgICAgIGlmKGlkeCA+PSBmcm9tX2NlbGxzLmxlbmd0aCB8fCBjb21wYXJlQ2VsbHMoZnJvbV9jZWxsc1tpZHhdLCBiKSAhPT0gMCkge1xuICAgICAgICAgIGJyZWFrXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgcmV0dXJuIGluZGV4XG59XG5leHBvcnRzLmluY2lkZW5jZSA9IGluY2lkZW5jZVxuXG4vL0NvbXB1dGVzIHRoZSBkdWFsIG9mIHRoZSBtZXNoLiAgVGhpcyBpcyBiYXNpY2FsbHkgYW4gb3B0aW1pemVkIHZlcnNpb24gb2YgYnVpbGRJbmRleCBmb3IgdGhlIHNpdHVhdGlvbiB3aGVyZSBmcm9tX2NlbGxzIGlzIGp1c3QgdGhlIGxpc3Qgb2YgdmVydGljZXNcbmZ1bmN0aW9uIGR1YWwoY2VsbHMsIHZlcnRleF9jb3VudCkge1xuICBpZighdmVydGV4X2NvdW50KSB7XG4gICAgcmV0dXJuIGluY2lkZW5jZSh1bmlxdWUoc2tlbGV0b24oY2VsbHMsIDApKSwgY2VsbHMsIDApXG4gIH1cbiAgdmFyIHJlcyA9IG5ldyBBcnJheSh2ZXJ0ZXhfY291bnQpXG4gIGZvcih2YXIgaT0wOyBpPHZlcnRleF9jb3VudDsgKytpKSB7XG4gICAgcmVzW2ldID0gW11cbiAgfVxuICBmb3IodmFyIGk9MCwgbGVuPWNlbGxzLmxlbmd0aDsgaTxsZW47ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MCwgY2w9Yy5sZW5ndGg7IGo8Y2w7ICsraikge1xuICAgICAgcmVzW2Nbal1dLnB1c2goaSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJlc1xufVxuZXhwb3J0cy5kdWFsID0gZHVhbFxuXG4vL0VudW1lcmF0ZXMgYWxsIGNlbGxzIGluIHRoZSBjb21wbGV4XG5mdW5jdGlvbiBleHBsb2RlKGNlbGxzKSB7XG4gIHZhciByZXN1bHQgPSBbXVxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgICAsIGNsID0gYy5sZW5ndGh8MFxuICAgIGZvcih2YXIgaj0xLCBqbD0oMTw8Y2wpOyBqPGpsOyArK2opIHtcbiAgICAgIHZhciBiID0gW11cbiAgICAgIGZvcih2YXIgaz0wOyBrPGNsOyArK2spIHtcbiAgICAgICAgaWYoKGogPj4+IGspICYgMSkge1xuICAgICAgICAgIGIucHVzaChjW2tdKVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXN1bHQucHVzaChiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gbm9ybWFsaXplKHJlc3VsdClcbn1cbmV4cG9ydHMuZXhwbG9kZSA9IGV4cGxvZGVcblxuLy9FbnVtZXJhdGVzIGFsbCBvZiB0aGUgbi1jZWxscyBvZiBhIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gc2tlbGV0b24oY2VsbHMsIG4pIHtcbiAgaWYobiA8IDApIHtcbiAgICByZXR1cm4gW11cbiAgfVxuICB2YXIgcmVzdWx0ID0gW11cbiAgICAsIGswICAgICA9ICgxPDwobisxKSktMVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGs9azA7IGs8KDE8PGMubGVuZ3RoKTsgaz1iaXRzLm5leHRDb21iaW5hdGlvbihrKSkge1xuICAgICAgdmFyIGIgPSBuZXcgQXJyYXkobisxKVxuICAgICAgICAsIGwgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajxjLmxlbmd0aDsgKytqKSB7XG4gICAgICAgIGlmKGsgJiAoMTw8aikpIHtcbiAgICAgICAgICBiW2wrK10gPSBjW2pdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJlc3VsdC5wdXNoKGIpXG4gICAgfVxuICB9XG4gIHJldHVybiBub3JtYWxpemUocmVzdWx0KVxufVxuZXhwb3J0cy5za2VsZXRvbiA9IHNrZWxldG9uO1xuXG4vL0NvbXB1dGVzIHRoZSBib3VuZGFyeSBvZiBhbGwgY2VsbHMsIGRvZXMgbm90IHJlbW92ZSBkdXBsaWNhdGVzXG5mdW5jdGlvbiBib3VuZGFyeShjZWxscykge1xuICB2YXIgcmVzID0gW11cbiAgZm9yKHZhciBpPTAsaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTAsY2w9Yy5sZW5ndGg7IGo8Y2w7ICsraikge1xuICAgICAgdmFyIGIgPSBuZXcgQXJyYXkoYy5sZW5ndGgtMSlcbiAgICAgIGZvcih2YXIgaz0wLCBsPTA7IGs8Y2w7ICsraykge1xuICAgICAgICBpZihrICE9PSBqKSB7XG4gICAgICAgICAgYltsKytdID0gY1trXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXMucHVzaChiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gbm9ybWFsaXplKHJlcylcbn1cbmV4cG9ydHMuYm91bmRhcnkgPSBib3VuZGFyeTtcblxuLy9Db21wdXRlcyBjb25uZWN0ZWQgY29tcG9uZW50cyBmb3IgYSBkZW5zZSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGNvbm5lY3RlZENvbXBvbmVudHNfZGVuc2UoY2VsbHMsIHZlcnRleF9jb3VudCkge1xuICB2YXIgbGFiZWxzID0gbmV3IFVuaW9uRmluZCh2ZXJ0ZXhfY291bnQpXG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wOyBqPGMubGVuZ3RoOyArK2opIHtcbiAgICAgIGZvcih2YXIgaz1qKzE7IGs8Yy5sZW5ndGg7ICsraykge1xuICAgICAgICBsYWJlbHMubGluayhjW2pdLCBjW2tdKVxuICAgICAgfVxuICAgIH1cbiAgfVxuICB2YXIgY29tcG9uZW50cyA9IFtdXG4gICAgLCBjb21wb25lbnRfbGFiZWxzID0gbGFiZWxzLnJhbmtzXG4gIGZvcih2YXIgaT0wOyBpPGNvbXBvbmVudF9sYWJlbHMubGVuZ3RoOyArK2kpIHtcbiAgICBjb21wb25lbnRfbGFiZWxzW2ldID0gLTFcbiAgfVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBsID0gbGFiZWxzLmZpbmQoY2VsbHNbaV1bMF0pXG4gICAgaWYoY29tcG9uZW50X2xhYmVsc1tsXSA8IDApIHtcbiAgICAgIGNvbXBvbmVudF9sYWJlbHNbbF0gPSBjb21wb25lbnRzLmxlbmd0aFxuICAgICAgY29tcG9uZW50cy5wdXNoKFtjZWxsc1tpXS5zbGljZSgwKV0pXG4gICAgfSBlbHNlIHtcbiAgICAgIGNvbXBvbmVudHNbY29tcG9uZW50X2xhYmVsc1tsXV0ucHVzaChjZWxsc1tpXS5zbGljZSgwKSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIGNvbXBvbmVudHNcbn1cblxuLy9Db21wdXRlcyBjb25uZWN0ZWQgY29tcG9uZW50cyBmb3IgYSBzcGFyc2UgZ3JhcGhcbmZ1bmN0aW9uIGNvbm5lY3RlZENvbXBvbmVudHNfc3BhcnNlKGNlbGxzKSB7XG4gIHZhciB2ZXJ0aWNlcyAgPSB1bmlxdWUobm9ybWFsaXplKHNrZWxldG9uKGNlbGxzLCAwKSkpXG4gICAgLCBsYWJlbHMgICAgPSBuZXcgVW5pb25GaW5kKHZlcnRpY2VzLmxlbmd0aClcbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTA7IGo8Yy5sZW5ndGg7ICsraikge1xuICAgICAgdmFyIHZqID0gZmluZENlbGwodmVydGljZXMsIFtjW2pdXSlcbiAgICAgIGZvcih2YXIgaz1qKzE7IGs8Yy5sZW5ndGg7ICsraykge1xuICAgICAgICBsYWJlbHMubGluayh2aiwgZmluZENlbGwodmVydGljZXMsIFtjW2tdXSkpXG4gICAgICB9XG4gICAgfVxuICB9XG4gIHZhciBjb21wb25lbnRzICAgICAgICA9IFtdXG4gICAgLCBjb21wb25lbnRfbGFiZWxzICA9IGxhYmVscy5yYW5rc1xuICBmb3IodmFyIGk9MDsgaTxjb21wb25lbnRfbGFiZWxzLmxlbmd0aDsgKytpKSB7XG4gICAgY29tcG9uZW50X2xhYmVsc1tpXSA9IC0xXG4gIH1cbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgbCA9IGxhYmVscy5maW5kKGZpbmRDZWxsKHZlcnRpY2VzLCBbY2VsbHNbaV1bMF1dKSk7XG4gICAgaWYoY29tcG9uZW50X2xhYmVsc1tsXSA8IDApIHtcbiAgICAgIGNvbXBvbmVudF9sYWJlbHNbbF0gPSBjb21wb25lbnRzLmxlbmd0aFxuICAgICAgY29tcG9uZW50cy5wdXNoKFtjZWxsc1tpXS5zbGljZSgwKV0pXG4gICAgfSBlbHNlIHtcbiAgICAgIGNvbXBvbmVudHNbY29tcG9uZW50X2xhYmVsc1tsXV0ucHVzaChjZWxsc1tpXS5zbGljZSgwKSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIGNvbXBvbmVudHNcbn1cblxuLy9Db21wdXRlcyBjb25uZWN0ZWQgY29tcG9uZW50cyBmb3IgYSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGNvbm5lY3RlZENvbXBvbmVudHMoY2VsbHMsIHZlcnRleF9jb3VudCkge1xuICBpZih2ZXJ0ZXhfY291bnQpIHtcbiAgICByZXR1cm4gY29ubmVjdGVkQ29tcG9uZW50c19kZW5zZShjZWxscywgdmVydGV4X2NvdW50KVxuICB9XG4gIHJldHVybiBjb25uZWN0ZWRDb21wb25lbnRzX3NwYXJzZShjZWxscylcbn1cbmV4cG9ydHMuY29ubmVjdGVkQ29tcG9uZW50cyA9IGNvbm5lY3RlZENvbXBvbmVudHNcbiIsIlwidXNlIHN0cmljdFwiXG5cbmZ1bmN0aW9uIHVuaXF1ZV9wcmVkKGxpc3QsIGNvbXBhcmUpIHtcbiAgdmFyIHB0ciA9IDFcbiAgICAsIGxlbiA9IGxpc3QubGVuZ3RoXG4gICAgLCBhPWxpc3RbMF0sIGI9bGlzdFswXVxuICBmb3IodmFyIGk9MTsgaTxsZW47ICsraSkge1xuICAgIGIgPSBhXG4gICAgYSA9IGxpc3RbaV1cbiAgICBpZihjb21wYXJlKGEsIGIpKSB7XG4gICAgICBpZihpID09PSBwdHIpIHtcbiAgICAgICAgcHRyKytcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGxpc3RbcHRyKytdID0gYVxuICAgIH1cbiAgfVxuICBsaXN0Lmxlbmd0aCA9IHB0clxuICByZXR1cm4gbGlzdFxufVxuXG5mdW5jdGlvbiB1bmlxdWVfZXEobGlzdCkge1xuICB2YXIgcHRyID0gMVxuICAgICwgbGVuID0gbGlzdC5sZW5ndGhcbiAgICAsIGE9bGlzdFswXSwgYiA9IGxpc3RbMF1cbiAgZm9yKHZhciBpPTE7IGk8bGVuOyArK2ksIGI9YSkge1xuICAgIGIgPSBhXG4gICAgYSA9IGxpc3RbaV1cbiAgICBpZihhICE9PSBiKSB7XG4gICAgICBpZihpID09PSBwdHIpIHtcbiAgICAgICAgcHRyKytcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGxpc3RbcHRyKytdID0gYVxuICAgIH1cbiAgfVxuICBsaXN0Lmxlbmd0aCA9IHB0clxuICByZXR1cm4gbGlzdFxufVxuXG5mdW5jdGlvbiB1bmlxdWUobGlzdCwgY29tcGFyZSwgc29ydGVkKSB7XG4gIGlmKGxpc3QubGVuZ3RoID09PSAwKSB7XG4gICAgcmV0dXJuIGxpc3RcbiAgfVxuICBpZihjb21wYXJlKSB7XG4gICAgaWYoIXNvcnRlZCkge1xuICAgICAgbGlzdC5zb3J0KGNvbXBhcmUpXG4gICAgfVxuICAgIHJldHVybiB1bmlxdWVfcHJlZChsaXN0LCBjb21wYXJlKVxuICB9XG4gIGlmKCFzb3J0ZWQpIHtcbiAgICBsaXN0LnNvcnQoKVxuICB9XG4gIHJldHVybiB1bmlxdWVfZXEobGlzdClcbn1cblxubW9kdWxlLmV4cG9ydHMgPSB1bmlxdWVcbiIsIlwidXNlIHN0cmljdFwiXG5cbnZhciBjaCA9IHJlcXVpcmUoXCJpbmNyZW1lbnRhbC1jb252ZXgtaHVsbFwiKVxudmFyIHVuaXEgPSByZXF1aXJlKFwidW5pcVwiKVxuXG5tb2R1bGUuZXhwb3J0cyA9IHRyaWFuZ3VsYXRlXG5cbmZ1bmN0aW9uIExpZnRlZFBvaW50KHAsIGkpIHtcbiAgdGhpcy5wb2ludCA9IHBcbiAgdGhpcy5pbmRleCA9IGlcbn1cblxuZnVuY3Rpb24gY29tcGFyZUxpZnRlZChhLCBiKSB7XG4gIHZhciBhcCA9IGEucG9pbnRcbiAgdmFyIGJwID0gYi5wb2ludFxuICB2YXIgZCA9IGFwLmxlbmd0aFxuICBmb3IodmFyIGk9MDsgaTxkOyArK2kpIHtcbiAgICB2YXIgcyA9IGJwW2ldIC0gYXBbaV1cbiAgICBpZihzKSB7XG4gICAgICByZXR1cm4gc1xuICAgIH1cbiAgfVxuICByZXR1cm4gMFxufVxuXG5mdW5jdGlvbiB0cmlhbmd1bGF0ZTFEKG4sIHBvaW50cywgaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICBpZihuID09PSAxKSB7XG4gICAgaWYoaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICAgICAgcmV0dXJuIFsgWy0xLCAwXSBdXG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBbXVxuICAgIH1cbiAgfVxuICB2YXIgbGlmdGVkID0gcG9pbnRzLm1hcChmdW5jdGlvbihwLCBpKSB7XG4gICAgcmV0dXJuIFsgcFswXSwgaSBdXG4gIH0pXG4gIGxpZnRlZC5zb3J0KGZ1bmN0aW9uKGEsYikge1xuICAgIHJldHVybiBhWzBdIC0gYlswXVxuICB9KVxuICB2YXIgY2VsbHMgPSBuZXcgQXJyYXkobiAtIDEpXG4gIGZvcih2YXIgaT0xOyBpPG47ICsraSkge1xuICAgIHZhciBhID0gbGlmdGVkW2ktMV1cbiAgICB2YXIgYiA9IGxpZnRlZFtpXVxuICAgIGNlbGxzW2ktMV0gPSBbIGFbMV0sIGJbMV0gXVxuICB9XG4gIGlmKGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgICBjZWxscy5wdXNoKFxuICAgICAgWyAtMSwgY2VsbHNbMF1bMV0sIF0sXG4gICAgICBbIGNlbGxzW24tMV1bMV0sIC0xIF0pXG4gIH1cbiAgcmV0dXJuIGNlbGxzXG59XG5cbmZ1bmN0aW9uIHRyaWFuZ3VsYXRlKHBvaW50cywgaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICB2YXIgbiA9IHBvaW50cy5sZW5ndGhcbiAgaWYobiA9PT0gMCkge1xuICAgIHJldHVybiBbXVxuICB9XG4gIFxuICB2YXIgZCA9IHBvaW50c1swXS5sZW5ndGhcbiAgaWYoZCA8IDEpIHtcbiAgICByZXR1cm4gW11cbiAgfVxuXG4gIC8vU3BlY2lhbCBjYXNlOiAgRm9yIDFEIHdlIGNhbiBqdXN0IHNvcnQgdGhlIHBvaW50c1xuICBpZihkID09PSAxKSB7XG4gICAgcmV0dXJuIHRyaWFuZ3VsYXRlMUQobiwgcG9pbnRzLCBpbmNsdWRlUG9pbnRBdEluZmluaXR5KVxuICB9XG4gIFxuICAvL0xpZnQgcG9pbnRzLCBzb3J0XG4gIHZhciBsaWZ0ZWQgPSBuZXcgQXJyYXkobilcbiAgdmFyIHVwcGVyID0gMS4wXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIHZhciBwID0gcG9pbnRzW2ldXG4gICAgdmFyIHggPSBuZXcgQXJyYXkoZCsxKVxuICAgIHZhciBsID0gMC4wXG4gICAgZm9yKHZhciBqPTA7IGo8ZDsgKytqKSB7XG4gICAgICB2YXIgdiA9IHBbal1cbiAgICAgIHhbal0gPSB2XG4gICAgICBsICs9IHYgKiB2XG4gICAgfVxuICAgIHhbZF0gPSBsXG4gICAgbGlmdGVkW2ldID0gbmV3IExpZnRlZFBvaW50KHgsIGkpXG4gICAgdXBwZXIgPSBNYXRoLm1heChsLCB1cHBlcilcbiAgfVxuICB1bmlxKGxpZnRlZCwgY29tcGFyZUxpZnRlZClcbiAgXG4gIC8vRG91YmxlIHBvaW50c1xuICBuID0gbGlmdGVkLmxlbmd0aFxuXG4gIC8vQ3JlYXRlIG5ldyBsaXN0IG9mIHBvaW50c1xuICB2YXIgZHBvaW50cyA9IG5ldyBBcnJheShuICsgZCArIDEpXG4gIHZhciBkaW5kZXggPSBuZXcgQXJyYXkobiArIGQgKyAxKVxuXG4gIC8vQWRkIHN0ZWluZXIgcG9pbnRzIGF0IHRvcFxuICB2YXIgdSA9IChkKzEpICogKGQrMSkgKiB1cHBlclxuICB2YXIgeSA9IG5ldyBBcnJheShkKzEpXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB5W2ldID0gMC4wXG4gIH1cbiAgeVtkXSA9IHVcblxuICBkcG9pbnRzWzBdID0geS5zbGljZSgpXG4gIGRpbmRleFswXSA9IC0xXG5cbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB4ID0geS5zbGljZSgpXG4gICAgeFtpXSA9IDFcbiAgICBkcG9pbnRzW2krMV0gPSB4XG4gICAgZGluZGV4W2krMV0gPSAtMVxuICB9XG5cbiAgLy9Db3B5IHJlc3Qgb2YgdGhlIHBvaW50cyBvdmVyXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIHZhciBoID0gbGlmdGVkW2ldXG4gICAgZHBvaW50c1tpICsgZCArIDFdID0gaC5wb2ludFxuICAgIGRpbmRleFtpICsgZCArIDFdID0gIGguaW5kZXhcbiAgfVxuXG4gIC8vQ29uc3RydWN0IGNvbnZleCBodWxsXG4gIHZhciBodWxsID0gY2goZHBvaW50cywgZmFsc2UpXG4gIGlmKGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgICBodWxsID0gaHVsbC5maWx0ZXIoZnVuY3Rpb24oY2VsbCkge1xuICAgICAgdmFyIGNvdW50ID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICB2YXIgdiA9IGRpbmRleFtjZWxsW2pdXVxuICAgICAgICBpZih2IDwgMCkge1xuICAgICAgICAgIGlmKCsrY291bnQgPj0gMikge1xuICAgICAgICAgICAgcmV0dXJuIGZhbHNlXG4gICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIGNlbGxbal0gPSB2XG4gICAgICB9XG4gICAgICByZXR1cm4gdHJ1ZVxuICAgIH0pXG4gIH0gZWxzZSB7XG4gICAgaHVsbCA9IGh1bGwuZmlsdGVyKGZ1bmN0aW9uKGNlbGwpIHtcbiAgICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgICAgdmFyIHYgPSBkaW5kZXhbY2VsbFtpXV1cbiAgICAgICAgaWYodiA8IDApIHtcbiAgICAgICAgICByZXR1cm4gZmFsc2VcbiAgICAgICAgfVxuICAgICAgICBjZWxsW2ldID0gdlxuICAgICAgfVxuICAgICAgcmV0dXJuIHRydWVcbiAgICB9KVxuICB9XG5cbiAgaWYoZCAmIDEpIHtcbiAgICBmb3IodmFyIGk9MDsgaTxodWxsLmxlbmd0aDsgKytpKSB7XG4gICAgICB2YXIgaCA9IGh1bGxbaV1cbiAgICAgIHZhciB4ID0gaFswXVxuICAgICAgaFswXSA9IGhbMV1cbiAgICAgIGhbMV0gPSB4XG4gICAgfVxuICB9XG5cbiAgcmV0dXJuIGh1bGxcbn0iLCJcbm1vZHVsZS5leHBvcnRzID0gcGFyc2VcblxuLyoqXG4gKiBleHBlY3RlZCBhcmd1bWVudCBsZW5ndGhzXG4gKiBAdHlwZSB7T2JqZWN0fVxuICovXG5cbnZhciBsZW5ndGggPSB7YTogNywgYzogNiwgaDogMSwgbDogMiwgbTogMiwgcTogNCwgczogNCwgdDogMiwgdjogMSwgejogMH1cblxuLyoqXG4gKiBzZWdtZW50IHBhdHRlcm5cbiAqIEB0eXBlIHtSZWdFeHB9XG4gKi9cblxudmFyIHNlZ21lbnQgPSAvKFthc3R2enFtaGxjXSkoW15hc3R2enFtaGxjXSopL2lnXG5cbi8qKlxuICogcGFyc2UgYW4gc3ZnIHBhdGggZGF0YSBzdHJpbmcuIEdlbmVyYXRlcyBhbiBBcnJheVxuICogb2YgY29tbWFuZHMgd2hlcmUgZWFjaCBjb21tYW5kIGlzIGFuIEFycmF5IG9mIHRoZVxuICogZm9ybSBgW2NvbW1hbmQsIGFyZzEsIGFyZzIsIC4uLl1gXG4gKlxuICogQHBhcmFtIHtTdHJpbmd9IHBhdGhcbiAqIEByZXR1cm4ge0FycmF5fVxuICovXG5cbmZ1bmN0aW9uIHBhcnNlKHBhdGgpIHtcblx0dmFyIGRhdGEgPSBbXVxuXHRwYXRoLnJlcGxhY2Uoc2VnbWVudCwgZnVuY3Rpb24oXywgY29tbWFuZCwgYXJncyl7XG5cdFx0dmFyIHR5cGUgPSBjb21tYW5kLnRvTG93ZXJDYXNlKClcblx0XHRhcmdzID0gcGFyc2VWYWx1ZXMoYXJncylcblxuXHRcdC8vIG92ZXJsb2FkZWQgbW92ZVRvXG5cdFx0aWYgKHR5cGUgPT0gJ20nICYmIGFyZ3MubGVuZ3RoID4gMikge1xuXHRcdFx0ZGF0YS5wdXNoKFtjb21tYW5kXS5jb25jYXQoYXJncy5zcGxpY2UoMCwgMikpKVxuXHRcdFx0dHlwZSA9ICdsJ1xuXHRcdFx0Y29tbWFuZCA9IGNvbW1hbmQgPT0gJ20nID8gJ2wnIDogJ0wnXG5cdFx0fVxuXG5cdFx0d2hpbGUgKHRydWUpIHtcblx0XHRcdGlmIChhcmdzLmxlbmd0aCA9PSBsZW5ndGhbdHlwZV0pIHtcblx0XHRcdFx0YXJncy51bnNoaWZ0KGNvbW1hbmQpXG5cdFx0XHRcdHJldHVybiBkYXRhLnB1c2goYXJncylcblx0XHRcdH1cblx0XHRcdGlmIChhcmdzLmxlbmd0aCA8IGxlbmd0aFt0eXBlXSkgdGhyb3cgbmV3IEVycm9yKCdtYWxmb3JtZWQgcGF0aCBkYXRhJylcblx0XHRcdGRhdGEucHVzaChbY29tbWFuZF0uY29uY2F0KGFyZ3Muc3BsaWNlKDAsIGxlbmd0aFt0eXBlXSkpKVxuXHRcdH1cblx0fSlcblx0cmV0dXJuIGRhdGFcbn1cblxuZnVuY3Rpb24gcGFyc2VWYWx1ZXMoYXJncyl7XG5cdGFyZ3MgPSBhcmdzLm1hdGNoKC8tP1suMC05XSsoPzplWy0rXT9cXGQrKT8vaWcpXG5cdHJldHVybiBhcmdzID8gYXJncy5tYXAoTnVtYmVyKSA6IFtdXG59XG4iLCIndXNlIHN0cmljdCc7XG5cbnZhciBzaWduID0gcmVxdWlyZSgnLi91dGlsaXRpZXMuanMnKS5zaWduO1xudmFyIGNhbGN1bGF0ZURpc3RhbmNlID0gcmVxdWlyZSgnLi91dGlsaXRpZXMuanMnKS5kaXN0YW5jZTtcblxuLy8gdmFyIHBvaW50cyA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnBvaW50cztcbi8vIHZhciBjaXR5U2V0ID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykuY2l0eVNldDtcbi8vIHZhciB0ZXh0UG9pbnRzSWQgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS50ZXh0UG9pbnRzSWQ7XG4vLyB2YXIgcG9zc2libGVTdGFydFBvaW50c0lkID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykucG9zc2libGVTdGFydFBvaW50c0lkO1xuXG52YXIgbGl2ZU1vdXNlUG9zaXRpb24gPSByZXF1aXJlKCcuL21vdXNlLmpzJyk7XG5cbnZhciBWZWN0b3IgPSByZXF1aXJlKCcuL3ZlY3Rvci5qcycpO1xuXG52YXIgcmFuZG9tID0gTWF0aC5yYW5kb207XG52YXIgZmxvb3IgPSBNYXRoLmZsb29yO1xuLy8gdmFyIFJFUFVMU0lPTiA9IDAuMDU7XG4vLyB2YXIgUkVQVUxTSU9OU1BFRUQgPSAwLjAwMjtcbi8vIHZhciBBTlRWRUxPQ0lUWSA9IDAuMDAxO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGNvbnRhaW5lciwgaW5pdFBvaW50cywgb3B0aW9ucyl7XG5cbiAgICBjb25zb2xlLmxvZygnT3B0aW9ucyBhbnQgOicsIG9wdGlvbnMpO1xuICAgIC8vIERlZmluZSB0aG9zZSBwYXJhbWV0ZXJzIGFzIGF0dHJpYnV0ZXMgb2YgQW50IG9iamVjdCA/XG4gICAgdmFyIFJFUFVMU0lPTiA9IG9wdGlvbnMucmVwU2l6ZTtcbiAgICB2YXIgUkVQVUxTSU9OU1BFRUQgPSBvcHRpb25zLnJlcFNwZWVkO1xuICAgIHZhciBBTlRWRUxPQ0lUWSA9IG9wdGlvbnMudmVsb2NpdHk7XG4gICAgdmFyIFdFSUdIVCA9IG9wdGlvbnMud2VpZ2h0O1xuXG4gICAgdmFyIG1vdXNlID0gbGl2ZU1vdXNlUG9zaXRpb24oY29udGFpbmVyKTtcblxuICAgIHZhciBwb2ludHMgPSBpbml0UG9pbnRzLnBvaW50cztcbiAgICB2YXIgY2l0eVNldCA9IGluaXRQb2ludHMuY2l0eVNldDtcbiAgICB2YXIgdGV4dFBvaW50c0lkID0gaW5pdFBvaW50cy50ZXh0UG9pbnRzSWQ7XG4gICAgdmFyIHBvc3NpYmxlU3RhcnRQb2ludHNJZCA9IGluaXRQb2ludHMucG9zc2libGVTdGFydFBvaW50c0lkO1xuXG5cbiAgICBmdW5jdGlvbiBBbnQocG9pbnQpIHtcbiAgICAgICAgdGhpcy54ID0gcG9pbnQueDsgICAgICAgICAgICAgICAgXG4gICAgICAgIHRoaXMueSA9IHBvaW50Lnk7XG4gICAgICAgIHRoaXMudmVsb2NpdHkgPSBBTlRWRUxPQ0lUWTtcbiAgICAgICAgdGhpcy53ZWlnaHQgPSBXRUlHSFQ7XG4gICAgICAgIHRoaXMucmVwU2l6ZSA9IFJFUFVMU0lPTjtcbiAgICAgICAgdGhpcy5yZXBTcGVlZCA9IFJFUFVMU0lPTlNQRUVEO1xuICAgICAgICB0aGlzLmVkZ2UgPSB1bmRlZmluZWQ7XG4gICAgICAgIHRoaXMuc3RhdGUgPSBcImZvcmFnZVwiO1xuICAgICAgICB0aGlzLmVkZ2VzID0gW107XG4gICAgICAgIHRoaXMubGFzdENpdHkgPSB1bmRlZmluZWQ7XG4gICAgICAgIHRoaXMub3JpZ2luID0gcG9pbnQ7XG4gICAgICAgIHRoaXMuZGVzdGluYXRpb24gPSB1bmRlZmluZWQ7XG4gICAgICAgIHRoaXMub3JpZW50YXRpb24gPSB1bmRlZmluZWQ7XG4gICAgICAgIHRoaXMuZGlyZWN0aW9uID0gbmV3IFZlY3RvcigwLDApO1xuICAgICAgICB0aGlzLnByb2cgPSAwO1xuICAgIH1cbiAgICAvLyBmb3JhZ2U6IHRoZSBhbnQgd2FuZGVycyBhcm91bmQgd2l0aG91dCBhbnkgcGhlcm9tb24gZGVwb3NpdGlvblxuICAgIC8vIG9uY2UgaXQgZmluZHMgYSBjaXR5LCBpdCBzdGFydHMgcmVtZW1iZXJpbmcgdGhlIG5vZGVzIGl0IGdvZXMgdGhyb3VnaFxuICAgIC8vIHdoZW4gaXQgZmluZHMgYW5vdGhlciBjaXR5LCBpdCBjb21wdXRlcyB0aGUgcGF0aCBsZW5ndGggYW5kIGFkZHMgcGhlcm9tb25zIG9uZSBlYWNoIGVkZ2VzXG4gICAgLy8gcHJvcG9ydGlvbm5hbHkgdG8gdGhlIHNob3J0ZXN0bmVzcyBvZiB0aGUgcGF0aFxuICAgIC8vIGl0IHJlc2V0cyB0aGUgbGlzdCBvZiBub2RlcyBhbmQgY29udGludWVzXG4gICAgLy8gd2hpbGUgZm9yYWdpbmcgdGhlIGFudCBjaG9zZXMgdGhlIHBhdGggd2l0aCBhIHBoZXJvbW9uIHByZWZlcmVuY2VcblxuXG4gICAgLy8gc3RhdGljIG1ldGhvZHNcbiAgICBBbnQuZ2VuZXJhdGVSYW5kU3RhcnRQb2ludCA9IGZ1bmN0aW9uKCkge1xuICAgICAgICB2YXIgcmFuZElkID0gTWF0aC5mbG9vcihwb3NzaWJsZVN0YXJ0UG9pbnRzSWQubGVuZ3RoICogcmFuZG9tKCkpO1xuICAgICAgICB2YXIgcmFuZFN0YXJ0UG9pbnQgPSBwb2ludHNbcG9zc2libGVTdGFydFBvaW50c0lkW3JhbmRJZF1dO1xuICAgICAgICByZXR1cm4gcmFuZFN0YXJ0UG9pbnQ7XG4gICAgfVxuXG5cbiAgICAvLyBtZXRob2RzXG4gICAgQW50LnByb3RvdHlwZSA9IHtcblxuICAgICAgICB0cmFuc2l0OiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgc3dpdGNoICh0aGlzLnN0YXRlKSB7XG4gICAgICAgICAgICBjYXNlIFwiZm9yYWdlXCI6XG4gICAgICAgICAgICAgICAgdmFyIHJlcyA9IHRoaXMubW92ZSgpO1xuICAgICAgICAgICAgICAgIGlmIChyZXMuY2l0eVJlYWNoZWQpIHtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5zdGF0ZSA9IFwicGhlcm9tb25pbmdcIjtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5sYXN0Q2l0eSA9IHRoaXMub3JpZ2luLmlkO1xuICAgICAgICAgICAgICAgIH07XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBjYXNlIFwicGhlcm9tb25pbmdcIjpcbiAgICAgICAgICAgICAgICB2YXIgcmVzID0gdGhpcy5tb3ZlKCk7XG4gICAgICAgICAgICAgICAgaWYgKHJlcy5lZGdlQ2hhbmdlZCkge1xuICAgICAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzLnB1c2godGhpcy5lZGdlKTtcbiAgICAgICAgICAgICAgICAgICAgLy8gZm91bmQgYSBjaXR5XG4gICAgICAgICAgICAgICAgICAgIGlmIChyZXMuY2l0eVJlYWNoZWQgJiYgKHRoaXMub3JpZ2luLmlkICE9IHRoaXMubGFzdENpdHkpICl7XG4gICAgICAgICAgICAgICAgICAgICAgICAvLyBjb21wdXRlIHRoZSBsZW5ndGggb2YgdGhlIHBhdGhcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBwYXRoTGVuZ3RoID0gdGhpcy5lZGdlcy5tYXAoZnVuY3Rpb24oZSl7cmV0dXJuIGUuZGlzdGFuY2V9KS5yZWR1Y2UoZnVuY3Rpb24oYSxiKXtyZXR1cm4gYSArIGJ9KTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBkZWx0YVBoZXJvbW9uZSA9IDEvcGF0aExlbmd0aDtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBhbnRXZWlnaHQgPSB0aGlzLndlaWdodDtcbiAgICAgICAgICAgICAgICAgICAgICAgIHRoaXMuZWRnZXMuZm9yRWFjaChmdW5jdGlvbihlKXtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgYSA9IGUucHQxLCBiID0gZS5wdDIsIHdlaWdodCA9IDE7ICBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAvLyBpbmNyZWFzZWQgZHJvcHBlZCBwaGVyb21vbnMgZm9yIHRleHRFZGdlc1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmICgoY2l0eVNldC5pbmRleE9mKGEuaWQpICE9IC0xKSAmJiBjaXR5U2V0LmluZGV4T2YoYi5pZCkgIT0gLTEgJiYgKE1hdGguYWJzKGEuaWQgLSBiLmlkKSA9PSAxKSlcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHdlaWdodCAqPSBhbnRXZWlnaHQ7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGUucGhlcm9tb24gKz0gKGRlbHRhUGhlcm9tb25lICogd2VpZ2h0KTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH0pO1xuXG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzID0gW3RoaXMuZWRnZV07XG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmxhc3RDaXR5ID0gdGhpcy5vcmlnaW4uaWQ7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgIH0sXG5cbiAgICAgICAgc2V0RGlyZWN0aW9uOiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgdmFyIHBvc3NpYmxlRWRnZXMgPSBbXTtcblxuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCB0aGlzLm9yaWdpbi5uZXh0cy5sZW5ndGg7IGkrKylcbiAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICBwb3NzaWJsZUVkZ2VzW2ldID0gdGhpcy5vcmlnaW4ubmV4dHNbaV07XG4gICAgICAgICAgICB9IFxuXG4gICAgICAgICAgICAvLyBjb25zb2xlLmxvZygnc21lbGxzMTogJywgcG9zc2libGVFZGdlcyk7XG5cbiAgICAgICAgICAgIHBvc3NpYmxlRWRnZXMuc3BsaWNlKHBvc3NpYmxlRWRnZXMuaW5kZXhPZih0aGlzLmVkZ2UpLDEpO1xuXG4gICAgICAgICAgICAvLyBmbGlwIGEgY29pbiBhbmQgZWl0aGVyIHRha2UgdGhlIHNtZWxsaWVzdCBwYXRoIG9yIGEgcmFuZG9tIG9uZVxuICAgICAgICAgICAgaWYgKHJhbmRvbSgpIDwgMC43NSl7XG4gICAgICAgICAgICAgICAgdmFyIHNtZWxscyA9IHBvc3NpYmxlRWRnZXMubWFwKGZ1bmN0aW9uKGUpe3JldHVybiBlLnBoZXJvbW9uO30pO1xuICAgICAgICAgICAgICAgIHZhciBpbmRleCA9IHNtZWxscy5pbmRleE9mKE1hdGgubWF4LmFwcGx5KE1hdGgsIHNtZWxscykpO1xuICAgICAgICAgICAgICAgIHRoaXMuZWRnZSA9IHBvc3NpYmxlRWRnZXNbaW5kZXhdO1xuICAgICAgICAgICAgfSBcbiAgICAgICAgICAgIGVsc2V7XG4gICAgICAgICAgICAgICAgdGhpcy5lZGdlID0gcG9zc2libGVFZGdlc1tmbG9vcihyYW5kb20oKSpwb3NzaWJsZUVkZ2VzLmxlbmd0aCldO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIFxuXG4gICAgICAgICAgICAvLyBzZXQgdGhlIGRlc3RpbmF0aW9uIHBvaW50LCBiZWluZyBlZGdlLnB0MSBvciBlZGdlLnB0MlxuICAgICAgICAgICAgdGhpcy5kZXN0aW5hdGlvbiA9ICh0aGlzLm9yaWdpbiA9PSB0aGlzLmVkZ2UucHQxKSA/IHRoaXMuZWRnZS5wdDIgOiB0aGlzLmVkZ2UucHQxO1xuXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi54ID0gdGhpcy5kZXN0aW5hdGlvbi54IC0gdGhpcy5vcmlnaW4ueDsgXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi55ID0gdGhpcy5kZXN0aW5hdGlvbi55IC0gdGhpcy5vcmlnaW4ueTtcblxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ubm9ybWFsaXplKCk7XG4gICAgICAgIH0sXG5cbiAgICAgICAgbW92ZTogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIC8vIGNvbnNvbGUubG9nKCdtb3ZlJyk7XG4gICAgICAgICAgICB2YXIgZWRnZUNoYW5nZWQ7XG4gICAgICAgICAgICB2YXIgY2l0eVJlYWNoZWQgPSBmYWxzZTtcblxuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ueCA9IHRoaXMuZGVzdGluYXRpb24ueCAtIHRoaXMueDsgXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi55ID0gdGhpcy5kZXN0aW5hdGlvbi55IC0gdGhpcy55O1xuICAgICAgICAgICAgdGhpcy5kaXJlY3Rpb24ubm9ybWFsaXplKCk7XG5cbiAgICAgICAgICAgIC8vIG9uIGVkZ2VcbiAgICAgICAgICAgIGlmICgoY2FsY3VsYXRlRGlzdGFuY2UodGhpcywgdGhpcy5kZXN0aW5hdGlvbikgPiB0aGlzLnJlcFNwZWVkKSl7XG5cbiAgICAgICAgICAgICAgICAvLyBhIGRlbHRhIG1vdmVtZW50IHdpbGwgYmUgYXBwbGllZCBpZiBjb2xsaXNpb24gd2l0aCBvYnN0YWNsZSBkZXRlY3RlZFxuICAgICAgICAgICAgICAgIHZhciBkZWx0YSA9IHRoaXMuYXZvaWRPYnN0YWNsZSgpO1xuXG4gICAgICAgICAgICAgICAgdGhpcy54ICs9IHRoaXMudmVsb2NpdHkgKiB0aGlzLmRpcmVjdGlvbi54ICsgZGVsdGEueCAqIHRoaXMucmVwU3BlZWQ7XG4gICAgICAgICAgICAgICAgdGhpcy55ICs9IHRoaXMudmVsb2NpdHkgKiB0aGlzLmRpcmVjdGlvbi55ICsgZGVsdGEueSAqIHRoaXMucmVwU3BlZWQ7XG5cbiAgICAgICAgICAgICAgICB0aGlzLnByb2cgPSB0aGlzLmNhbGN1bGF0ZVByb2dyZXNzaW9uKCk7XG4gICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgZWRnZUNoYW5nZWQgPSBmYWxzZTtcblxuICAgICAgICAgICAgLy8gb24gdmVydGV4XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIC8vIGNvbnNvbGUubG9nKCdyZWFjaGVkJyk7XG4gICAgICAgICAgICAgICAgdGhpcy5zdGVwID0gMDtcbiAgICAgICAgICAgICAgICB0aGlzLnByb2cgPSAwO1xuICAgICAgICAgICAgICAgIHRoaXMub3JpZ2luID0gdGhpcy5kZXN0aW5hdGlvbjtcbiAgICAgICAgICAgICAgICB0aGlzLnggPSB0aGlzLm9yaWdpbi54O1xuICAgICAgICAgICAgICAgIHRoaXMueSA9IHRoaXMub3JpZ2luLnk7XG5cbiAgICAgICAgICAgICAgICB0aGlzLnNldERpcmVjdGlvbigpO1xuXG4gICAgICAgICAgICAgICAgY2l0eVJlYWNoZWQgPSAoY2l0eVNldC5pbmRleE9mKHRoaXMub3JpZ2luLmlkKSAhPSAtMSk7XG4gICAgICAgICAgICAgICAgZWRnZUNoYW5nZWQgPSB0cnVlO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmV0dXJuIHtjaXR5UmVhY2hlZDogY2l0eVJlYWNoZWQsIGVkZ2VDaGFuZ2VkOiBlZGdlQ2hhbmdlZH07XG4gICAgICAgIH0sXG5cbiAgICAgICAgYXZvaWRPYnN0YWNsZTogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHZhciBkaXN0YW5jZSA9IGNhbGN1bGF0ZURpc3RhbmNlKHRoaXMsIG1vdXNlKTtcbiAgICAgICAgXG4gICAgICAgICAgICBpZiAoZGlzdGFuY2UgPD0gdGhpcy5yZXBTaXplKSB7XG5cbiAgICAgICAgICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgICAgICAgICAvLyBkZWx0YSBtb3ZlbWVudCBpcyBjb21wb3NlZCBvZiBhIHJlcHVsc2lvbiBkZWx0YSBhbmQgYSBjaXJjdWxhciBkZWx0YSBcbiAgICAgICAgICAgICAgICAgICAgeDogKHRoaXMueCAtIG1vdXNlLngpL2Rpc3RhbmNlICsgKHRoaXMueSAtIG1vdXNlLnkpL2Rpc3RhbmNlICogMSxcbiAgICAgICAgICAgICAgICAgICAgeTogKHRoaXMueSAtIG1vdXNlLnkpL2Rpc3RhbmNlIC0gKHRoaXMueCAtIG1vdXNlLngpL2Rpc3RhbmNlICogMVxuICAgICAgICAgICAgICAgIH07XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBlbHNlXG4gICAgICAgICAgICAgICAgcmV0dXJuIHt4OjAsIHk6MH07XG4gICAgICAgIH0sXG5cbiAgICAgICAgY2FsY3VsYXRlUHJvZ3Jlc3Npb246IGZ1bmN0aW9uKCl7XG4gICAgICAgICAgICB2YXIgdiA9IG5ldyBWZWN0b3IodGhpcy54IC0gdGhpcy5vcmlnaW4ueCwgdGhpcy55IC0gdGhpcy5vcmlnaW4ueSk7XG4gICAgICAgICAgICB2YXIgbm9ybSA9IHYubm9ybSgpO1xuXG4gICAgICAgICAgICB2YXIgdGhldGEgPSAodi54ICogdGhpcy5lZGdlLmRpcmVjdGlvbi54ICsgdi55ICogdGhpcy5lZGdlLmRpcmVjdGlvbi55KSAvIG5vcm07XG4gICAgICAgICAgICB2YXIgcHJvZyA9IG5vcm0gKiBNYXRoLmFicyh0aGV0YSk7XG4gICAgICAgICAgICAvLyByZXR1cm5zIGxlbmd0aCBvZiBwcm9qZWN0aW9uIG9uIGVkZ2VcbiAgICAgICAgICAgIHJldHVybiBwcm9nO1xuICAgICAgICB9XG5cbiAgICB9O1xuICAgIHJldHVybiBBbnQ7XG59XG5cbiIsIid1c2Ugc3RyaWN0J1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uIChBbnQpIHtcblxuXHR2YXIgbmJBbnRzUGVyU3RlcCA9IDEwMDtcblxuXHRmdW5jdGlvbiBjcmVhdGVHcm91cChwb3B1bGF0aW9uKXtcblx0XHRmb3IgKHZhciBpID0gMDsgaSA8IG5iQW50c1BlclN0ZXA7IGkrKykge1xuXHRcdFx0dmFyIG5ld0FudCA9IG5ldyBBbnQoQW50LmdlbmVyYXRlUmFuZFN0YXJ0UG9pbnQoKSk7XG5cdFx0XHRuZXdBbnQuc2V0RGlyZWN0aW9uKCk7XG5cdFx0XHRwb3B1bGF0aW9uLnB1c2gobmV3QW50KTtcblx0XHR9XG5cbi8vIFx0XHRjb25zb2xlLmxvZygnQ3JlYXRlZCBBbnRzIEdyb3VwOiBcXFxuLy8gKCsgJyArIG5iQW50c1BlclN0ZXAgKyAnKSA9PiAnICsgcG9wdWxhdGlvbi5sZW5ndGgpO1xuXG5cdFx0cmV0dXJuIHBvcHVsYXRpb247XG5cdH1cblxuXHRmdW5jdGlvbiByZW1vdmVHcm91cChwb3B1bGF0aW9uLCBuYkRlYWQpe1xuXHRcdHBvcHVsYXRpb24gPSBwb3B1bGF0aW9uLnNsaWNlKDAsIHBvcHVsYXRpb24ubGVuZ3RoIC0gbmJEZWFkKTtcblxuLy8gXHRcdGNvbnNvbGUubG9nKCdSZW1vdmVkIEFudHMgR3JvdXA6IFxcXG4vLyAoLSAnICsgbmJBbnRzUGVyU3RlcCArICcpID0+ICcgKyBwb3B1bGF0aW9uLmxlbmd0aCk7XG5cblx0XHRyZXR1cm4gcG9wdWxhdGlvbjtcblxuXHR9XG5cblx0cmV0dXJuIHtcblx0XHRjcmVhdGU6IGNyZWF0ZUdyb3VwLFxuXHRcdHJlbW92ZTogcmVtb3ZlR3JvdXBcblx0fTtcblxufVxuXHQiLCIndXNlIHN0cmljdCdcblxudmFyIGR0ID0gcmVxdWlyZShcImRlbGF1bmF5LXRyaWFuZ3VsYXRlXCIpO1xuXG52YXIgcmFuZ2UgPSByZXF1aXJlKCcuL3V0aWxpdGllcy5qcycpLnJhbmdlO1xuXG4vLyB2YXIgcG9pbnRzID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykucG9pbnRzO1xudmFyIHRleHRNZXNoID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykudGV4dE1lc2g7XG52YXIgY2l0eVNldCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLmNpdHlTZXQ7XG52YXIgbmJSYW5kb21Qb2ludHMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5uYlJhbmRvbVBvaW50cztcbnZhciBmb3JjZWRFZGdlcyA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLmZvcmNlZEVkZ2VzO1xuXG52YXIgRWRnZSA9IHJlcXVpcmUoJy4vZWRnZS5qcycpO1xuXG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24ocG9pbnRzKXtcbiAgICAvLyB0cmlhbmd1bGF0ZVxuICAgIHZhciBjZWxscyA9IGR0KHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7XG4gICAgICAgIHJldHVybiBbcC54LCBwLnldO1xuICAgIH0pKTtcblxuICAgIHZhciBlZGdlcyA9IFtdO1xuICAgIHZhciBwZXJtdXRhdGlvbnMgPSBbWzAsMV0sIFswLDJdLCBbMSwyXV07XG5cbiAgICAvLyBmb3JjZSB0aGUgZWRnZXMgb2YgdGhlIHRleHQgdG8gYmUgZWRnZXMgb2YgdGhlIGdyYXBoXG4gICAgaWYgKHRleHRNZXNoKSB7XG4gICAgICAgIHJhbmdlKDAsIHBvaW50cy5sZW5ndGggLSBuYlJhbmRvbVBvaW50cykuZm9yRWFjaChmdW5jdGlvbihpZCl7XG4gICAgICAgICAgICB2YXIgZGlyZWN0TGluayA9IGZvcmNlZEVkZ2VzW2lkXTtcbiAgICAgICAgICAgIHZhciB0ZXh0RWRnZSA9IEVkZ2UuY3JlYXRlKHBvaW50c1tpZF0sIHBvaW50c1tkaXJlY3RMaW5rXSk7XG4gICAgICAgICAgICBlZGdlcy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgIHBvaW50c1tpZF0ubmV4dHMucHVzaCh0ZXh0RWRnZSk7XG4gICAgICAgIH0pXG4gICAgfVxuXG5cbiAgICBjZWxscy5mb3JFYWNoKGZ1bmN0aW9uKGNlbGwpe1xuICAgICAgIFxuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IDM7ICsraSl7ICAvLyBmb3IgZWFjaCBwb2ludC5pZCBsaXN0ZWQgaW4gY3VycmVudCBjZWxsXG4gICAgICAgICAgICB2YXIgcHQgPSBwb2ludHNbY2VsbFtpXV07XG5cbiAgICAgICAgICAgIGZvciAodmFyIGogPSAxOyBqIDwgMzsgKytqKXsgXG5cbiAgICAgICAgICAgICAgICB2YXIgcHRqID0gcG9pbnRzW2NlbGxbKCBpICsgaiApICUgM11dOyAvLyBwaWNrIG9uZSBvZiB0aGUgb3RoZXIgMiBwb2ludHMgb2YgdGhlIGNlbGxcbiAgICAgICAgICAgICAgICB2YXIgbmV3RWRnZSA9IHVuZGVmaW5lZDtcblxuICAgICAgICAgICAgICAgIC8vIGlmIHB0IGFscmVhZHkgaGFzIG5leHRFZGdlcyAuLi5cbiAgICAgICAgICAgICAgICBpZiAocHQubmV4dHMubGVuZ3RoICE9IDApIHtcbiAgICAgICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgICAgIC8vIC4uLiBnZXQgdGhlIHBvaW50cyBjb3JyZXNwb25kaW5nIC4uLlxuICAgICAgICAgICAgICAgICAgICB2YXIgdGVtcFBvaW50cyA9IHB0Lm5leHRzLm1hcChmdW5jdGlvbihlKXtcbiAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBbZS5wdDEsIGUucHQyXTtcbiAgICAgICAgICAgICAgICAgICAgfSkucmVkdWNlKGZ1bmN0aW9uKGEsIGIpe1xuICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBhLmNvbmNhdChiKTtcbiAgICAgICAgICAgICAgICAgICAgfSk7XG5cbiAgICAgICAgICAgICAgICAgICAgLy8gLi4uIGFuZCBjaGVjayBpZiBwdGogYWxyZWFkeSBpcyBwYXJ0IG9mIHRoZSBleGlzdGluZyBuZXh0RWRnZXMuIElmIG5vdCwgYWRkIHRoZSBlZGdlLlxuICAgICAgICAgICAgICAgICAgICBpZiAodGVtcFBvaW50cy5pbmRleE9mKHB0aikgPT0gLTEpe1xuICAgICAgICAgICAgICAgICAgICAgICAgbmV3RWRnZSA9IEVkZ2UuY3JlYXRlKHB0LCBwdGopO1xuICAgICAgICAgICAgICAgICAgICAgICAgZWRnZXMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHB0Lm5leHRzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIG5ld0VkZ2UgPSBFZGdlLmNyZWF0ZShwdCwgcHRqKTtcbiAgICAgICAgICAgICAgICAgICAgZWRnZXMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICAgICAgcHQubmV4dHMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAvLyBhZGQgYWxzbyB0aGUgZWRnZSB0byB0aGUgZWRnZSdzIG90aGVyIHBvaW50J3MgbmV4dEVkZ2VzXG4gICAgICAgICAgICAgICAgaWYgKG5ld0VkZ2UgIT0gdW5kZWZpbmVkKXtcbiAgICAgICAgICAgICAgICAgICAgcHRqLm5leHRzLnB1c2gobmV3RWRnZSk7XG4gICAgICAgICAgICAgICAgfSAgICAgICAgIFxuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAvLyBhZGQgdGhlIHRleHRFZGdlcyB0byBuZXh0RWRnZXMgbWFwXG4gICAgICAgICAgICBpZiAodGV4dE1lc2ggJiYgKGNpdHlTZXQuaW5kZXhPZihwdCkgIT0gLTEpKSB7XG4gICAgICAgICAgICAgICAgdmFyIHRleHRFZGdlID0gRWRnZS5jcmVhdGUocHQsIHBvaW50c1twdC5pZCArIDFdKTtcbiAgICAgICAgICAgICAgICBlZGdlcy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgICAgICBwdC5uZXh0cy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICB9XG4gICAgfSk7XG5cbiAgICByZXR1cm4gZWRnZXM7XG59OyIsIid1c2Ugc3RyaWN0JztcblxudmFyIHNxcnQgPSBNYXRoLnNxcnQ7XG52YXIgcG93ID0gTWF0aC5wb3c7XG52YXIgYWJzID0gTWF0aC5hYnM7XG52YXIgYXRhbiA9IE1hdGguYXRhbjtcblxudmFyIFZlY3RvciA9IHJlcXVpcmUoJy4vdmVjdG9yLmpzJyk7XG5cblxuZnVuY3Rpb24gRWRnZShwdEEsIHB0Qikge1xuICAgIHZhciBkaXN0YW5jZSA9IHNxcnQoIHBvdyhwdEEueCAtIHB0Qi54LCAyKSArIHBvdyhwdEEueSAtIHB0Qi55LCAyKSApO1xuXG4gICAgLy8gZmluZCBsaW5lIGVxdWF0aW9uIGF4ICsgYnkgKyBjID0gMFxuICAgIHZhciBhID0gMTtcbiAgICB2YXIgYiA9IC0gKHB0Qi54IC0gcHRBLngpIC8gKHB0Qi55IC0gcHRBLnkpO1xuXG4gICAgLy8gb3JpZW50YXRlIHZlY3RvciAoYSxiKVxuICAgIGlmIChiIDwgMCl7XG4gICAgICAgIGIgPSAtYjtcbiAgICAgICAgYSA9IC1hO1xuICAgIH1cblxuICAgIC8vIG5vcm1hbGl6ZSB2ZWN0b3IgKGEsYilcbiAgICB2YXIgbiA9IG5ldyBWZWN0b3IoYSwgYik7XG4gICAgbi5ub3JtYWxpemUoKTtcblxuICAgIHZhciBjID0gLSAoYSAqIHB0QS54ICsgYiAqIHB0QS55KTtcblxuICAgIC8vIC8vIGNhbGN1bGF0ZSB2ZWN0b3IgZGlyZWN0b3JcbiAgICB2YXIgdiA9IG5ldyBWZWN0b3IocHRCLnggLSBwdEEueCwgcHRCLnkgLSBwdEEueSk7XG4gICAgXG4gICAgdi5ub3JtYWxpemUoKTtcblxuICAgIHRoaXMuaWQgPSB1bmRlZmluZWQ7XG4gICAgdGhpcy5wdDEgPSBwdEE7XG4gICAgdGhpcy5wdDIgPSBwdEI7XG4gICAgdGhpcy5kaXJlY3Rpb24gPSB2O1xuICAgIHRoaXMub3J0aERpcmVjdGlvbiA9IG47IFxuICAgIHRoaXMuZGlzdGFuY2UgPSBkaXN0YW5jZTtcbiAgICB0aGlzLnBoZXJvbW9uID0gMS9kaXN0YW5jZTtcbiAgICB0aGlzLmxpbmUgPSB7XG4gICAgICAgIGE6IGEsXG4gICAgICAgIGI6IGIsXG4gICAgICAgIGM6IGMsXG4gICAgfTtcblxuICAgIGlmICh0aGlzLmRpc3RhbmNlID09PSAwKSBjb25zb2xlLmxvZygnWkVSTyAhJyk7XG59XG5cblxuLy8gc3RhdGljIG1ldGhvZHNcbkVkZ2UuY3JlYXRlID0gZnVuY3Rpb24ocHRBLCBwdEIpIHtcbiAgICB2YXIgZWRnZSA9IG5ldyBFZGdlKHB0QSwgcHRCKTtcbiAgICByZXR1cm4gZWRnZTtcbn1cblxuXG4vLyBtZXRob2RzXG5FZGdlLnByb3RvdHlwZSA9IHtcblxuICAgIGdldE90aGVyUG9pbnQ6IGZ1bmN0aW9uKHBvaW50KSB7XG4gICAgICAgIGlmIChwb2ludCA9PSB0aGlzLnB0MSlcbiAgICAgICAgICAgIHJldHVybiB0aGlzLnB0MjtcbiAgICAgICAgZWxzZSBpZiAocG9pbnQgPT0gdGhpcy5wdDIpXG4gICAgICAgICAgICByZXR1cm4gdGhpcy5wdDE7XG4gICAgICAgIGVsc2VcbiAgICAgICAgICAgIGNvbnNvbGUubG9nKFwiRXJyb3JcIik7XG4gICAgfSxcblxuICAgIGNhbGN1bGF0ZURpc3RhbmNlOiBmdW5jdGlvbih4LCB5KSB7XG4gICAgICAgIHZhciBhID0gdGhpcy5saW5lLmEsXG4gICAgICAgICAgICBiID0gdGhpcy5saW5lLmIsXG4gICAgICAgICAgICBjID0gdGhpcy5saW5lLmM7XG4gICAgICAgIHJldHVybiBhYnMoYSAqIHggKyBiICogeSArIGMpIC8gTWF0aC5zcXJ0KE1hdGgucG93KGEsMikgKyBNYXRoLnBvdyhiLDIpKTtcbiAgICB9LFxuXG59XG5tb2R1bGUuZXhwb3J0cyA9IEVkZ2U7IiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgcGFyc2UgPSByZXF1aXJlKCdwYXJzZS1zdmctcGF0aCcpO1xuXG52YXIgcmFuZ2UgPSByZXF1aXJlKCcuL3V0aWxpdGllcy5qcycpLnJhbmdlO1xuXG52YXIgUG9pbnQgPSByZXF1aXJlKCcuL3BvaW50LmpzJyk7XG52YXIgc3ZnUGF0aCA9IHJlcXVpcmUoJy4vc3ZnUGF0aC5qcycpO1xuXG52YXIgcmFuZG9tID0gTWF0aC5yYW5kb207XG5cbnZhciBuYkNpdHkgPSAyO1xuXG52YXIgdGV4dE1lc2ggPSB0cnVlO1xuXG4vLyBGcmFtZSBkZWZpbml0aW9uXG52YXIgeEluaXQgPSAwLCB5SW5pdCA9IDA7XG52YXIgdyA9IDEsXG4gICAgaCA9IDE7XG5cbnZhciBzdmdTdHJpbmcgPSBzdmdQYXRoO1xuXG5mdW5jdGlvbiBzdmdUb1BvaW50cyhzdmdTdHJpbmcpIHtcbiAgICB2YXIgcG9pbnRzID0gW107XG4gICAgdmFyIGVkZ2VzID0gT2JqZWN0LmNyZWF0ZShudWxsKTtcblxuICAgIHZhciBiZWdpbmluZ1BhdGg7XG5cbiAgICB2YXIgWCA9IDA7XG4gICAgdmFyIFkgPSAwO1xuICAgIHZhciBuYlBvaW50cyA9IDA7XG4gICAgdmFyIHByZXZQb2ludDtcblxuICAgIHZhciBjb21tYW5kcyA9IHBhcnNlKHN2Z1N0cmluZylcbiAgICBmb3IgKHZhciBpPTA7IGk8Y29tbWFuZHMubGVuZ3RoOyBpKyspe1xuICAgICAgICB2YXIgY29tbWFuZCA9IGNvbW1hbmRzW2ldO1xuICAgICAgICBzd2l0Y2ggKGNvbW1hbmRbMF0pIHtcbiAgICAgICAgICAgIGNhc2UgXCJtXCI6XG4gICAgICAgICAgICAgICAgWCArPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgKz0gY29tbWFuZFsyXTtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYmVnaW5pbmdQYXRoID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICBjYXNlIFwiTVwiOlxuICAgICAgICAgICAgICAgIFggPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgPSBjb21tYW5kWzJdO1xuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBiZWdpbmluZ1BhdGggPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICBicmVhazsgXG4gICAgICAgICAgICBjYXNlIFwiY1wiOlxuICAgICAgICAgICAgICAgIFggPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgPSBjb21tYW5kWzJdO1xuICAgICAgICAgICAgICAgIHBvaW50cy5wdXNoKHtpZDpuYlBvaW50cywgeDpYLCB5Oll9KTtcblxuICAgICAgICAgICAgICAgIGlmIChwcmV2UG9pbnQgIT0gdW5kZWZpbmVkKSB7XG4gICAgICAgICAgICAgICAgICAgIGVkZ2VzW3ByZXZQb2ludF0gPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgcHJldlBvaW50ID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgbmJQb2ludHMrKztcbiAgICAgICAgICAgICAgICBicmVhazsgXG4gICAgICAgICAgICBjYXNlIFwibFwiOlxuICAgICAgICAgICAgICAgIFggPSBjb21tYW5kWzFdO1xuICAgICAgICAgICAgICAgIFkgPSBjb21tYW5kWzJdO1xuICAgICAgICAgICAgICAgIHBvaW50cy5wdXNoKHtpZDpuYlBvaW50cywgeDpYLCB5Oll9KTtcblxuICAgICAgICAgICAgICAgIGlmIChwcmV2UG9pbnQgIT0gdW5kZWZpbmVkKSB7XG4gICAgICAgICAgICAgICAgICAgIGVkZ2VzW3ByZXZQb2ludF0gPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgcHJldlBvaW50ID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgbmJQb2ludHMrKztcbiAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgIGNhc2UgXCJ6XCI6XG4gICAgICAgICAgICAgICAgZWRnZXNbcHJldlBvaW50XSA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIGJlZ2luaW5nUGF0aCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSB1bmRlZmluZWQ7XG4gICAgICAgICAgICAgICAgYnJlYWs7ICAgIFxuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiB7cG9pbnRzIDogcG9pbnRzLCBlZGdlcyA6IGVkZ2VzfTtcbn1cblxuLy8gaW5pdGlhbGl6ZSBwb2ludHNcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihuYlN0YXJ0UG9pbnRzLCBuYlJhbmRvbVBvaW50cyl7XG4gICAgdmFyIHBvaW50cyA9IFtdO1xuICAgIHZhciBmb3JjZWRFZGdlcztcbiAgICB2YXIgY2l0eVNldDtcblxuICAgIGlmICh0ZXh0TWVzaCl7XG5cbiAgICAgICAgdmFyIG15VGV4dCA9IHN2Z1RvUG9pbnRzKHN2Z1N0cmluZyk7XG4gICAgICAgIHBvaW50cyA9IG15VGV4dC5wb2ludHM7XG4gICAgICAgIGZvcmNlZEVkZ2VzID0gbXlUZXh0LmVkZ2VzO1xuICAgICAgICBjaXR5U2V0ID0gcmFuZ2UoMCwgcG9pbnRzLmxlbmd0aCk7XG5cbiAgICAgICAgdmFyIHNjYWxlWCA9IDAuNTtcbiAgICAgICAgdmFyIHNjYWxlWSA9IDAuNDtcbiAgICAgICAgdmFyIGRlbHRhWCA9IDAuMjU7XG4gICAgICAgIHZhciBkZWx0YVkgPSAwLjI1O1xuXG4gICAgICAgIC8vIHNjYWxlIHBvaW50cyB0byBbMCwxXSArIGRlbHRhXG4gICAgICAgIHZhciBtYXhYID0gTWF0aC5tYXguYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC54fSkpO1xuICAgICAgICB2YXIgbWluWCA9IE1hdGgubWluLmFwcGx5KE1hdGgsIHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIHAueH0pKTtcbiAgICAgICAgdmFyIG1heFkgPSBNYXRoLm1heC5hcHBseShNYXRoLCBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiBwLnl9KSk7XG4gICAgICAgIHZhciBtaW5ZID0gTWF0aC5taW4uYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC55fSkpO1xuICAgICAgICBwb2ludHMgPSBwb2ludHMubWFwKGZ1bmN0aW9uKHApe1xuICAgICAgICAgICAgdmFyIHggPSBzY2FsZVggKiAocC54LW1pblgpLyhtYXhYLW1pblgpICsgZGVsdGFYO1xuICAgICAgICAgICAgdmFyIHkgPSBzY2FsZVkgKiAocC55LW1pblkpLyhtYXhZLW1pblkpICsgZGVsdGFZO1xuICAgICAgICAgICAgdmFyIG5ld1BvaW50ID0gbmV3IFBvaW50KHgsIHkpO1xuICAgICAgICAgICAgbmV3UG9pbnQuaWQgPSBwLmlkO1xuXG4gICAgICAgICAgICByZXR1cm4gbmV3UG9pbnQ7XG4gICAgICAgIH0pO1xuXG4gICAgICAgIC8vIG9ubHkgYWRkIHJhbmRvbSBwb2ludHNcbiAgICAgICAgdmFyIG5iUG9pbnRzID0gcG9pbnRzLmxlbmd0aDtcbiAgICAgICAgZm9yKHZhciBpPTA7IGk8bmJSYW5kb21Qb2ludHM7ICsraSkge1xuXG4gICAgICAgICAgICB2YXIgeCA9IHJhbmRvbSgpO1xuICAgICAgICAgICAgdmFyIHkgPSByYW5kb20oKTtcblxuICAgICAgICAgICAgdmFyIG5ld1BvaW50ID0gbmV3IFBvaW50KHgsIHkpO1xuICAgICAgICAgICAgbmV3UG9pbnQuaWQgPSBuYlBvaW50cztcblxuICAgICAgICAgICAgcG9pbnRzLnB1c2gobmV3UG9pbnQpO1xuXG4gICAgICAgICAgICBuYlBvaW50cysrO1xuICAgICAgICB9XG5cbiAgICB9IGVsc2Uge1xuICAgICAgICAvL2FkZCByYW5kb20gcG9pbnRzXG5cbiAgICAgICAgdmFyIG5iUG9pbnRzID0gMDtcbiAgICAgICAgZm9yKHZhciBpPTA7IGk8bmJSYW5kb21Qb2ludHM7ICsraSkge1xuXG4gICAgICAgICAgICB2YXIgeCA9IHJhbmRvbSgpO1xuICAgICAgICAgICAgdmFyIHkgPSByYW5kb20oKTtcblxuICAgICAgICAgICAgdmFyIG5ld1BvaW50ID0gbmV3IFBvaW50KHgsIHkpO1xuICAgICAgICAgICAgbmV3UG9pbnQuaWQgPSBuYlBvaW50cztcblxuICAgICAgICAgICAgcG9pbnRzLnB1c2gobmV3UG9pbnQpO1xuICAgICAgICAgICAgXG4gICAgICAgICAgICBuYlBvaW50cysrO1xuICAgICAgICB9XG5cbiAgICAgICAgY2l0eVNldCA9IHJhbmdlKDAsIG5iQ2l0eSk7XG4gICAgICAgIGNvbnNvbGUubG9nKGNpdHlTZXQpO1xuICAgIH1cblxuXG4gICAgLy8gaW5pdGlhbGl6ZSBzdGFydCBwb2ludHNcbiAgICB2YXIgcG9zc2libGVTdGFydFBvaW50c0lkID0gW107XG5cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IG5iU3RhcnRQb2ludHM7IGkrKyl7XG4gICAgICAgIHBvc3NpYmxlU3RhcnRQb2ludHNJZC5wdXNoKE1hdGguZmxvb3IobmJQb2ludHMgKiByYW5kb20oKSkpO1xuICAgIH1cbiAgICBcblxuICAgIHJldHVybiB7XG4gICAgICAgIHRleHRNZXNoOiB0ZXh0TWVzaCxcbiAgICAgICAgcG9pbnRzOiBwb2ludHMsXG4gICAgICAgIGNpdHlTZXQ6IGNpdHlTZXQsXG4gICAgICAgIHBvc3NpYmxlU3RhcnRQb2ludHNJZDogcG9zc2libGVTdGFydFBvaW50c0lkLFxuICAgICAgICBuYlJhbmRvbVBvaW50czogbmJSYW5kb21Qb2ludHMsXG4gICAgICAgIGZvcmNlZEVkZ2VzOiBmb3JjZWRFZGdlc1xuICAgIH07XG59XG4iLCIndXNlIHN0cmljdCdcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbiAoY29udGFpbmVyKXtcblxuXHR2YXIgbW91c2UgPSB7XG5cdCAgICB4OiAwLFxuXHQgICAgeTogMFxuXHR9O1xuXG5cdGNvbnRhaW5lci5hZGRFdmVudExpc3RlbmVyKCAnbW91c2Vtb3ZlJywgZnVuY3Rpb24oZSl7XG5cdCAgICB2YXIgcmVjdCA9IGNvbnRhaW5lci5nZXRCb3VuZGluZ0NsaWVudFJlY3QoKTtcblxuXHQgICAgbW91c2UueCA9IChlLmNsaWVudFggLSByZWN0LmxlZnQgKSAvIHJlY3Qud2lkdGg7XG5cdCAgICBtb3VzZS55ID0gKGUuY2xpZW50WSAtIHJlY3QudG9wICkvIHJlY3QuaGVpZ2h0O1xuXHR9KTtcblxuXHRyZXR1cm4gbW91c2U7XG5cbn07XG4iLCIndXNlIHN0cmljdCdcblxuZnVuY3Rpb24gUG9pbnQoeCwgeSkge1xuICAgIHRoaXMuaWQgPSB1bmRlZmluZWQ7ICAgICAgICAgICAgICAgIFxuICAgIHRoaXMueCA9IHg7XG4gICAgdGhpcy55ID0geTtcbiAgICB0aGlzLm5leHRzID0gW107XG59XG5cbm1vZHVsZS5leHBvcnRzID0gUG9pbnQ7IiwiJ3VzZSBzdHJpY3QnXG5cbnZhciBhbnRGdW5jdGlvbiA9IHJlcXVpcmUoJy4vYW50LmpzJyk7XG52YXIgYW50c0dyb3VwID0gcmVxdWlyZSgnLi9hbnRzR3JvdXAnKTtcblxudmFyIHJhbmRvbSA9IE1hdGgucmFuZG9tO1xuXG52YXIgUkFORE9NTVZUID0gMC4wMDM7XG52YXIgQU5UU0laRSA9IDAuMDAyO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKGNvbnRhaW5lciwgcG9pbnRzTWFwLCBvcHRpb25zKXtcblxuXHRpZighY29udGFpbmVyKVxuXHRcdHRocm93IG5ldyBUeXBlRXJyb3IoJ01pc3NpbmcgY29udGFpbmVyJyk7XG5cblx0Ly8gQW50cyB2YXJpYWJsZXNcblx0dmFyIGVkZ2VzID0gcG9pbnRzTWFwLmVkZ2VzO1xuXHR2YXIgb2JqUG9wdWxhdGlvbkluaXRpYWwgPSBvcHRpb25zLm5iQW50cztcblx0dmFyIG9ialBvcHVsYXRpb24gPSBvYmpQb3B1bGF0aW9uSW5pdGlhbDtcblx0dmFyIHBvaW50c0luZm9zID0gcG9pbnRzTWFwLnBvaW50c0luZm9zO1xuXHR2YXIgcG9wdWxhdGlvbiA9IFtdO1xuXHR2YXIgbmJBbnRzUGVyU3RlcCA9IDEwMDtcblx0XG5cdHZhciBBbnQgPSBhbnRGdW5jdGlvbihjb250YWluZXIsIHBvaW50c0luZm9zLCBvcHRpb25zKTtcblx0YW50c0dyb3VwID0gYW50c0dyb3VwKEFudCk7XG5cblx0Ly8gQW5pbWF0aW9uIHZhcmlhYmxlc1xuXHR2YXIgYW5pbUlEO1xuXHR2YXIgZGVsdGFUaW1lO1xuXHR2YXIgRlBTQ291bnQ7XG5cdHZhciBsYXN0VXBkYXRlID0gcGVyZm9ybWFuY2Uubm93KCk7XG5cdHZhciBGUFNNb25pdG9yID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcignI0ZQUycpO1xuXHR2YXIgZFRNb25pdG9yID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcignI2RUJyk7XG5cdHZhciB3YXJuaW5nTW9uaXRvciA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3IoJyN3YXJuaW5nJyk7XG5cdHZhciByZWZyZXNoVGltZSA9IDA7XG5cdHZhciBtYXhEZWx0YVRpbWUgPSA0MDtcblx0dmFyIEZQU092ZXJMaW1pdENvdW50ID0gMDtcblx0dmFyIEZQU1VuZGVyTGltaXRDb3VudCA9IDA7XG5cblx0Ly8gd2FybmluZyBtZXNzYWdlIGRpc2FwcGVhcnMgYWZ0ZXIgNCBzXG5cdHdhcm5pbmdNb25pdG9yLmFkZEV2ZW50TGlzdGVuZXIoXCJ0cmFuc2l0aW9uZW5kXCIsIGZ1bmN0aW9uKCl7XG5cdFx0d2luZG93LnNldFRpbWVvdXQoZnVuY3Rpb24oKXtcblx0XHRcdHdhcm5pbmdNb25pdG9yLmNsYXNzTmFtZSA9IFwiaW52aXNpYmxlXCI7XG5cdFx0fSwgNDAwMCk7XG5cdH0pXG5cblx0Ly8gQ2FudmFzXG5cdHZhciBjYW52YXNMaXN0ID0gZG9jdW1lbnQuZ2V0RWxlbWVudHNCeVRhZ05hbWUoXCJjYW52YXNcIik7XG5cdFxuXHRpZiAoY2FudmFzTGlzdC5sZW5ndGggPT09IDApe1xuXHRcdHZhciBjYW52YXMgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50KFwiY2FudmFzXCIpO1xuXHRcdHZhciByZWN0ID0gY29udGFpbmVyLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpO1xuXHRcdGNhbnZhcy53aWR0aCA9IHJlY3Qud2lkdGg7XG5cdFx0Y2FudmFzLmhlaWdodCA9IHJlY3QuaGVpZ2h0O1xuXHRcdGNhbnZhcy5zdHlsZS5iYWNrZ3JvdW5kQ29sb3IgPSBcInJnYmEoMjUwLCAyNTAsIDI1MCwgMClcIjsgXG5cdFx0Y29udGFpbmVyLmFwcGVuZENoaWxkKGNhbnZhcyk7XG5cdH1cblx0ZWxzZXtcblx0XHR2YXIgY2FudmFzID0gY2FudmFzTGlzdFswXTtcblx0XHRjb25zb2xlLmxvZygnQ0FOVkFTJyk7XG5cdH1cblx0dmFyIGNvbnRleHQgPSBjYW52YXMuZ2V0Q29udGV4dChcIjJkXCIpO1xuXHRjb250ZXh0LmNsZWFyUmVjdCAoIDAgLCAwICwgY2FudmFzLndpZHRoLCBjYW52YXMuaGVpZ2h0ICk7XG5cdFxuXG5cdGZ1bmN0aW9uIGNoZWNrQW50TnVtYmVyKGFudE51bWJlcil7XG5cdFx0aWYgKGFudE51bWJlciA8IG9ialBvcHVsYXRpb24gLSA1MCl7XG5cdFx0XHRGUFNNb25pdG9yLnN0eWxlLmNvbG9yID0gXCJncmVlblwiO1xuXHRcdFx0cG9wdWxhdGlvbiA9IGFudHNHcm91cC5jcmVhdGUocG9wdWxhdGlvbik7XG5cdFx0fVx0XG5cdFx0ZWxzZSBpZiAoYW50TnVtYmVyID4gb2JqUG9wdWxhdGlvbil7XG5cdFx0XHRwb3B1bGF0aW9uID0gYW50c0dyb3VwLnJlbW92ZShwb3B1bGF0aW9uLCBhbnROdW1iZXIgLSBvYmpQb3B1bGF0aW9uKTtcblx0XHRcdEZQU01vbml0b3Iuc3R5bGUuY29sb3IgPSBcInJlZFwiO1xuXHRcdFx0d2FybmluZ01vbml0b3IuY2xhc3NOYW1lID0gXCJ2aXNpYmxlXCI7XG5cdFx0fVxuXHRcdGVsc2V7XG5cdFx0XHRGUFNNb25pdG9yLnN0eWxlLmNvbG9yID0gXCJ3aGl0ZVwiO1xuXHRcdFx0Ly8gd2FybmluZ01vbml0b3IuY2xhc3NOYW1lID0gXCJpbnZpc2libGVcIjtcblx0XHR9XG5cdH1cblxuXHRmdW5jdGlvbiBkaXNwbGF5RlBTKGRUKXtcblx0XHRGUFNDb3VudCA9ICgxMDAwL2RUKS50b0ZpeGVkKDIpO1xuXHRcdHZhciB0ID0gZFQudG9GaXhlZCgyKTtcblx0XHRGUFNNb25pdG9yLnRleHRDb250ZW50ID0gJ0ZQUyA6ICcgKyBGUFNDb3VudDsgIFxuXHRcdGRUTW9uaXRvci50ZXh0Q29udGVudCA9ICduYkFudHMgOiAnICsgcG9wdWxhdGlvbi5sZW5ndGg7XG5cdFx0Ly8gZFRNb25pdG9yLmlubmVyVGV4dCA9ICdkVCA6ICcgKyB0ICsgJ21zJztcblx0fVxuXG5cdGZ1bmN0aW9uIHRpY2soKSB7XG5cdFx0dmFyIG5vdyA9IHBlcmZvcm1hbmNlLm5vdygpO1xuXHRcdGRlbHRhVGltZSA9IG5vdyAtIGxhc3RVcGRhdGU7XG5cdFx0bGFzdFVwZGF0ZSA9IG5vdztcblx0XHRyZWZyZXNoVGltZSArPSBkZWx0YVRpbWUvMTAwMDsgLy8gaW4gc2Vjb25kc1xuXG5cdFx0Ly8gY29uc29sZS5sb2coJ25iQW50cycsIHBvcHVsYXRpb24ubGVuZ3RoKTtcblxuXHRcdGNoZWNrQW50TnVtYmVyKHBvcHVsYXRpb24ubGVuZ3RoKTtcblxuXHRcdC8vIGRpc3BsYXkgRlBTIGluZm8gZXZlcnkgMC4zIHNcblx0XHRpZiAocmVmcmVzaFRpbWUgPiAwLjMpe1xuXHRcdFx0ZGlzcGxheUZQUyhkZWx0YVRpbWUpO1xuXHRcdFx0cmVmcmVzaFRpbWUgPSAwOyBcblx0XHR9XG5cblx0XHQvLyByZW1vdmUgYW50cyB3aGVuIGZyYW1lIHJhdGUgaXMgdG9vIGxvd1xuXHRcdGlmIChGUFNPdmVyTGltaXRDb3VudCA9PT0gMTApIHtcblx0XHRcdG9ialBvcHVsYXRpb24gPSBvYmpQb3B1bGF0aW9uICogbWF4RGVsdGFUaW1lIC8gZGVsdGFUaW1lO1xuXHRcdFx0RlBTT3ZlckxpbWl0Q291bnQgPSAwO1xuXHRcdH1cblxuXHRcdHdoaWxlIChGUFNVbmRlckxpbWl0Q291bnQgPiA1MCAmJiBvYmpQb3B1bGF0aW9uIDwgb2JqUG9wdWxhdGlvbkluaXRpYWwpIHtcblx0XHRcdG9ialBvcHVsYXRpb24gKz0gMTA7XG5cdFx0fVxuXG5cdFx0Ly8gY2hlY2sgZHVyYXRpb24gb2Ygb3Zlci91bmRlciBmcmFtZXJhdGUgbGltaXQgcGVyaW9kc1xuXHRcdGlmIChkZWx0YVRpbWUgPiBtYXhEZWx0YVRpbWUpe1xuXHRcdFx0RlBTT3ZlckxpbWl0Q291bnQrKztcblx0XHRcdEZQU1VuZGVyTGltaXRDb3VudCA9IDA7XG5cdFx0fVxuXHRcdGVsc2Uge1xuXHRcdFx0RlBTT3ZlckxpbWl0Q291bnQgPSAwO1xuXHRcdFx0RlBTVW5kZXJMaW1pdENvdW50Kys7XG5cdFx0fVxuXG5cdFx0Ly8gZHJhdyBpbiBjYW52YXNcblx0XHR2YXIgdyA9IGNhbnZhcy53aWR0aDtcblx0XHR2YXIgaCA9IGNhbnZhcy5oZWlnaHQ7XG5cdFx0dmFyIG1vdXNlID0gW2xhc3RNb3VzZU1vdmVFdmVudC5jbGllbnRYL3csIGxhc3RNb3VzZU1vdmVFdmVudC5jbGllbnRZL2hdO1xuXHRcdGNvbnRleHQuc2V0VHJhbnNmb3JtKHcsIDAsIDAsIGgsIDAsIDApO1xuXHRcdGNvbnRleHQuZmlsbFN0eWxlID0gXCJyZ2JhKDI1MCwgMjUwLCAyNTAsIDAuNClcIjtcblx0XHRjb250ZXh0LmZpbGxSZWN0KDAsMCx3LGgpO1xuXG5cdFx0Ly8gLy8gZWRnZXNcblx0XHQvLyBjb250ZXh0LnN0cm9rZVN0eWxlID0gXCIjMDAwXCI7XG5cdFx0Ly8gZm9yKHZhciBpPTA7IGkgPCBlZGdlcy5sZW5ndGg7ICsraSkge1xuXHRcdC8vICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IDAuMDAwMTtcblx0XHQvLyAgICAgdmFyIGVkZ2UgPSBlZGdlc1tpXTtcblx0XHQvLyAgICAgaWYgKGVkZ2UucGhlcm9tb24gIT0gMCl7XG5cdFx0Ly8gICAgICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IE1hdGgubWluKDAuMDAwMDEgKiBlZGdlLnBoZXJvbW9uLCAwLjAxKTtcblx0XHQvLyAgICAgfSBlbHNlIHtcblx0XHQvLyAgICAgICAgIGNvbnRleHQubGluZVdpZHRoID0gMC4wMDAwMTtcblx0XHQvLyAgICAgfVxuXHRcdC8vICAgICBjb250ZXh0LmJlZ2luUGF0aCgpO1xuXHRcdC8vICAgICBjb250ZXh0Lm1vdmVUbyhwb2ludHNJbmZvcy5wb2ludHNbZWRnZS5wdDEuaWRdLngsIHBvaW50c0luZm9zLnBvaW50c1tlZGdlLnB0MS5pZF0ueSk7XG5cdFx0Ly8gICAgIGNvbnRleHQubGluZVRvKHBvaW50c0luZm9zLnBvaW50c1tlZGdlLnB0Mi5pZF0ueCwgcG9pbnRzSW5mb3MucG9pbnRzW2VkZ2UucHQyLmlkXS55KTtcblx0XHQvLyAgICAgY29udGV4dC5zdHJva2UoKTtcblx0XHQvLyB9XG5cblx0XHQvLyAvLyB2ZXJ0aWNlc1xuXHRcdC8vIGZvcih2YXIgaT0wOyBpPHBvaW50c0luZm9zLnBvaW50cy5sZW5ndGg7ICsraSkge1xuXHRcdC8vICAgICBjb250ZXh0LmJlZ2luUGF0aCgpXG5cdFx0Ly8gICAgIHZhciBwb2ludCA9IHBvaW50c0luZm9zLnBvaW50c1tpXTtcblx0XHQvLyAgICAgaWYgKHBvaW50c0luZm9zLmNpdHlTZXQuaW5kZXhPZihwb2ludC5pZCkgIT0gLTEpe1xuXHRcdC8vICAgICAgICAgY29udGV4dC5maWxsU3R5bGUgPSBcIiMwMTAxREZcIjtcblx0XHQvLyAgICAgICAgIGNvbnRleHQuYXJjKHBvaW50LngsIHBvaW50LnksIDAuMDA2LCAwLCAyKk1hdGguUEkpO1xuXHRcdC8vICAgICB9XG5cdFx0Ly8gICAgIGVsc2Uge1xuXHRcdC8vICAgICAgICAgY29udGV4dC5maWxsU3R5bGUgPSBcIiMwMDBcIjtcblx0XHQvLyAgICAgICAgIGNvbnRleHQuYXJjKHBvaW50c0luZm9zLnBvaW50c1tpXS54LCBwb2ludHNJbmZvcy5wb2ludHNbaV0ueSwgMC4wMDMsIDAsIDIqTWF0aC5QSSk7XG5cdFx0Ly8gICAgIH1cblx0XHQvLyAgICAgY29udGV4dC5jbG9zZVBhdGgoKTtcblx0XHQvLyAgICAgY29udGV4dC5maWxsKCk7XG5cdFx0Ly8gfVxuXG5cdFx0Ly8gbW92ZSBhbnRzXG5cdFx0cG9wdWxhdGlvbi5mb3JFYWNoKGZ1bmN0aW9uKGFudCl7XG5cdFx0XHRhbnQudHJhbnNpdCgpO1xuXHRcdH0pO1xuXG5cdFx0Ly8gcGhlcm9tb24gZXZhcG9yYXRpb25cblx0XHRlZGdlcy5mb3JFYWNoKGZ1bmN0aW9uKGVkZ2Upe1xuXHRcdFx0aWYoZWRnZS5waGVyb21vbiA+IDApe1xuXHRcdFx0XHRlZGdlLnBoZXJvbW9uIC09IDAuMDAwMTtcblx0XHRcdH1cblx0XHR9KTtcblxuXHRcdC8vIGFudHNcblx0XHRwb3B1bGF0aW9uLmZvckVhY2goZnVuY3Rpb24oYW50KXtcblx0XHRcdGNvbnRleHQuYmVnaW5QYXRoKClcblx0XHRcdHZhciB4ID0gYW50LnggKyBSQU5ET01NVlQqcmFuZG9tKCk7XG5cdFx0XHR2YXIgeSA9IGFudC55ICsgUkFORE9NTVZUKnJhbmRvbSgpO1xuXG5cdFx0XHRjb250ZXh0LmZpbGxTdHlsZSA9IFwiYmxhY2tcIlxuXHRcdFx0Y29udGV4dC5maWxsUmVjdCh4LCB5LCBBTlRTSVpFLCBBTlRTSVpFKTtcblx0XHRcdGNvbnRleHQuY2xvc2VQYXRoKCk7XG5cdFx0XHRjb250ZXh0LmZpbGwoKTtcblx0XHR9KVxuXHR9O1xuXHRcblx0dmFyIGxhc3RNb3VzZU1vdmVFdmVudCA9IHtcblx0XHRjbGllbnRYOiAwLFxuXHRcdGNsaWVudFk6IDBcblx0fTtcblx0XG5cdGNvbnRhaW5lci5hZGRFdmVudExpc3RlbmVyKCdtb3VzZW1vdmUnLCBmdW5jdGlvbihlKXtcblx0XHRsYXN0TW91c2VNb3ZlRXZlbnQgPSBlO1xuXHR9KTtcblx0XG5cdHZhciBwYXVzZWQgPSBmYWxzZTtcblx0XG5cdGZ1bmN0aW9uIHRvZ2dsZVBsYXlQYXVzZSgpe1xuXHRcdHBhdXNlZCA9ICFwYXVzZWQ7XG5cdFx0aWYoIXBhdXNlZClcblx0XHRcdGFuaW1hdGUoKTtcblx0fVxuXG5cdGZ1bmN0aW9uIHJlc2V0KCl7XG5cdFx0cG9wdWxhdGlvbiA9IFtdO1xuXHRcdGVkZ2VzID0gW107XG5cdFx0cG9pbnRzSW5mb3MgPSBbXTtcblxuXHRcdGNhbmNlbEFuaW1hdGlvbkZyYW1lKGFuaW1JRCk7XG5cdH1cblx0XG5cdC8vIGNvbnRhaW5lci5hZGRFdmVudExpc3RlbmVyKCdjbGljaycsIHRvZ2dsZVBsYXlQYXVzZSk7XG5cblx0ZnVuY3Rpb24gYW5pbWF0ZSgpe1xuXHRcdHRpY2soKTtcblx0XHRcblx0XHRpZighcGF1c2VkKVxuXHRcdFx0YW5pbUlEID0gcmVxdWVzdEFuaW1hdGlvbkZyYW1lKGFuaW1hdGUpO1xuXHR9XG5cdGFuaW1hdGUoKTtcblxuXHRmdW5jdGlvbiBtb2RpZnlBbnRzKG9wdHMpe1xuXHRcdG9ialBvcHVsYXRpb24gPSBvcHRzLm5iQW50cztcblxuXHRcdHBvcHVsYXRpb24uZm9yRWFjaChmdW5jdGlvbihhbnQpe1xuXHRcdFx0YW50LnZlbG9jaXR5ID0gb3B0cy52ZWxvY2l0eTtcblx0XHRcdGFudC53ZWlnaHQgPSBvcHRzLndlaWdodDtcblx0XHRcdGFudC5yZXBTaXplID0gb3B0cy5yZXBTaXplO1xuXHRcdFx0YW50LnJlcFNwZWVkID0gb3B0cy5yZXBTcGVlZDtcblx0XHR9KTtcblx0fVxuXHRcblx0cmV0dXJuIHtcblx0XHR0b2dnbGVQbGF5UGF1c2U6IHRvZ2dsZVBsYXlQYXVzZSxcblx0XHRyZXNldDogcmVzZXQsXG5cdFx0Ly8gc2hvdWxkIGJlIGEgZ2V0dGVyL3NldHRlciwgYnV0IElFOFxuXHRcdGdldEFudENvdW50OiBmdW5jdGlvbigpe1xuXHRcdFx0cmV0dXJuIHBvcHVsYXRpb24ubGVuZ3RoO1xuXHRcdH0sXG5cdFx0bW9kaWZ5QW50czogbW9kaWZ5QW50c1xuXHR9XG59XG4iLCIndXNlIHN0cmljdCc7XG5cbnZhciBwb2ludDEgPSBcIm0gMTguMjUsMTkuNSBjIDE4LjI1LDE1Ljc2Nzc2ODggMTUuMDAwMDk1MSwxMi43NSAxMSwxMi43NSBjIDYuOTk5OTA0ODgsMTIuNzUgMy43NSwxNS43Njc3Njg4IDMuNzUsMTkuNSBjIDMuNzUsMjMuMjMyMjMxMiA2Ljk5OTkwNDg4LDI2LjI1IDExLDI2LjI1IGMgMTUuMDAwMDk1MSwyNi4yNSAxOC4yNSwyMy4yMzIyMzEyIDE4LjI1LDE5LjUgeiBtIDQuMjUsMTkuNSBjIDQuMjUsMTYuMDUyNTI5NCA3LjI2ODEwODYzLDEzLjI1IDExLDEzLjI1IGMgMTQuNzMxODkxNCwxMy4yNSAxNy43NSwxNi4wNTI1Mjk0IDE3Ljc1LDE5LjUgYyAxNy43NSwyMi45NDc0NzA2IDE0LjczMTg5MTQsMjUuNzUgMTEsMjUuNzUgYyA3LjI2ODEwODYzLDI1Ljc1IDQuMjUsMjIuOTQ3NDcwNiA0LjI1LDE5LjUgelwiO1xudmFyIHBvaW50MiA9IFwibSA4OS4yNSw4LjUgYyA4OS4yNSw0LjIxNjA1Nzg3IDg1LjU1Mjg3MTYsMC43NSA4MSwwLjc1IGMgNzYuNDQ3MTI4NCwwLjc1IDcyLjc1LDQuMjE2MDU3ODcgNzIuNzUsOC41IGMgNzIuNzUsMTIuNzgzOTQyMSA3Ni40NDcxMjg0LDE2LjI1IDgxLDE2LjI1IGMgODUuNTUyODcxNiwxNi4yNSA4OS4yNSwxMi43ODM5NDIxIDg5LjI1LDguNSB6IG0gNzMuMjUsOC41IGMgNzMuMjUsNC40OTk2NzA4OCA3Ni43MTYzMTU2LDEuMjUgODEsMS4yNSBjIDg1LjI4MzY4NDQsMS4yNSA4OC43NSw0LjQ5OTY3MDg4IDg4Ljc1LDguNSBjIDg4Ljc1LDEyLjUwMDMyOTEgODUuMjgzNjg0NCwxNS43NSA4MSwxNS43NSBjIDc2LjcxNjMxNTYsMTUuNzUgNzMuMjUsMTIuNTAwMzI5MSA3My4yNSw4LjUgelwiO1xudmFyIHBvaW50MyA9IFwibSAxNjAuMjUsMTEgYyAxNjAuMjUsNS4zMzM0ODM3NSAxNTUuMjA4MTY4LDAuNzUgMTQ5LDAuNzUgYyAxNDIuNzkxODMyLDAuNzUgMTM3Ljc1LDUuMzMzNDgzNzUgMTM3Ljc1LDExIGMgMTM3Ljc1LDE2LjY2NjUxNjIgMTQyLjc5MTgzMiwyMS4yNSAxNDksMjEuMjUgYyAxNTUuMjA4MTY4LDIxLjI1IDE2MC4yNSwxNi42NjY1MTYyIDE2MC4yNSwxMSB6IG0gMTM4LjI1LDExIGMgMTM4LjI1LDUuNjIwODIxMjUgMTQzLjA1NzkwMywxLjI1IDE0OSwxLjI1IGMgMTU0Ljk0MjA5NywxLjI1IDE1OS43NSw1LjYyMDgyMTI1IDE1OS43NSwxMSBjIDE1OS43NSwxNi4zNzkxNzg3IDE1NC45NDIwOTcsMjAuNzUgMTQ5LDIwLjc1IGMgMTQzLjA1NzkwMywyMC43NSAxMzguMjUsMTYuMzc5MTc4NyAxMzguMjUsMTEgelwiO1xudmFyIHBvaW50NCA9IFwibSAxNjAuMjUsNzYuNSBjIDE2MC4yNSw3Mi43Njc3Njg4IDE1Ny4wMDAwOTUsNjkuNzUgMTUzLDY5Ljc1IGMgMTQ4Ljk5OTkwNSw2OS43NSAxNDUuNzUsNzIuNzY3NzY4OCAxNDUuNzUsNzYuNSBjIDE0NS43NSw4MC4yMzIyMzEyIDE0OC45OTk5MDUsODMuMjUgMTUzLDgzLjI1IGMgMTU3LjAwMDA5NSw4My4yNSAxNjAuMjUsODAuMjMyMjMxMiAxNjAuMjUsNzYuNSB6IG0gMTQ2LjI1LDc2LjUgYyAxNDYuMjUsNzMuMDUyNTI5NCAxNDkuMjY4MTA5LDcwLjI1IDE1Myw3MC4yNSBjIDE1Ni43MzE4OTEsNzAuMjUgMTU5Ljc1LDczLjA1MjUyOTQgMTU5Ljc1LDc2LjUgYyAxNTkuNzUsNzkuOTQ3NDcwNiAxNTYuNzMxODkxLDgyLjc1IDE1Myw4Mi43NSBjIDE0OS4yNjgxMDksODIuNzUgMTQ2LjI1LDc5Ljk0NzQ3MDYgMTQ2LjI1LDc2LjUgelwiO1xudmFyIHBvaW50NSA9IFwibSA5NS4yNSw3NiBjIDk1LjI1LDcwLjMzNjI3OTUgOTAuNDM0NDA2NSw2NS43NSA4NC41LDY1Ljc1IGMgNzguNTY1NTkzNSw2NS43NSA3My43NSw3MC4zMzYyNzk1IDczLjc1LDc2IGMgNzMuNzUsODEuNjYzNzIwNSA3OC41NjU1OTM1LDg2LjI1IDg0LjUsODYuMjUgYyA5MC40MzQ0MDY1LDg2LjI1IDk1LjI1LDgxLjY2MzcyMDUgOTUuMjUsNzYgeiBtIDc0LjI1LDc2IGMgNzQuMjUsNzAuNjE4MDI1NSA3OC44MzY0MjY4LDY2LjI1IDg0LjUsNjYuMjUgYyA5MC4xNjM1NzMyLDY2LjI1IDk0Ljc1LDcwLjYxODAyNTUgOTQuNzUsNzYgYyA5NC43NSw4MS4zODE5NzQ1IDkwLjE2MzU3MzIsODUuNzUgODQuNSw4NS43NSBjIDc4LjgzNjQyNjgsODUuNzUgNzQuMjUsODEuMzgxOTc0NSA3NC4yNSw3NiB6XCI7XG52YXIgcG9pbnQ2ID0gXCJtIDIwLjI1LDc1IGMgMjAuMjUsNzAuOTkxOTMzOCAxNi43NzY0OTk1LDY3Ljc1IDEyLjUsNjcuNzUgYyA4LjIyMzUwMDQ2LDY3Ljc1IDQuNzUsNzAuOTkxOTMzOCA0Ljc1LDc1IGMgNC43NSw3OS4wMDgwNjYyIDguMjIzNTAwNDYsODIuMjUgMTIuNSw4Mi4yNSBjIDE2Ljc3NjQ5OTUsODIuMjUgMjAuMjUsNzkuMDA4MDY2MiAyMC4yNSw3NSB6IG0gNS4yNSw3NSBjIDUuMjUsNzEuMjc2MDc5NyA4LjQ5MjIyODI5LDY4LjI1IDEyLjUsNjguMjUgYyAxNi41MDc3NzE3LDY4LjI1IDE5Ljc1LDcxLjI3NjA3OTcgMTkuNzUsNzUgYyAxOS43NSw3OC43MjM5MjAzIDE2LjUwNzc3MTcsODEuNzUgMTIuNSw4MS43NSBjIDguNDkyMjI4MjksODEuNzUgNS4yNSw3OC43MjM5MjAzIDUuMjUsNzUgelwiO1xudmFyIHBvaW50NyA9IFwibSAyMy4yNSwxMzkgYyAyMy4yNSwxMzIuNzg2Nzk3IDE4LjIxMzIwMzQsMTI3Ljc1IDEyLDEyNy43NSBjIDUuNzg2Nzk2NTYsMTI3Ljc1IDAuNzUsMTMyLjc4Njc5NyAwLjc1LDEzOSBjIDAuNzUsMTQ1LjIxMzIwMyA1Ljc4Njc5NjU2LDE1MC4yNSAxMiwxNTAuMjUgYyAxOC4yMTMyMDM0LDE1MC4yNSAyMy4yNSwxNDUuMjEzMjAzIDIzLjI1LDEzOSB6IG0gMS4yNSwxMzkgYyAxLjI1LDEzMy4wNjI5MzkgNi4wNjI5Mzg5NCwxMjguMjUgMTIsMTI4LjI1IGMgMTcuOTM3MDYxMSwxMjguMjUgMjIuNzUsMTMzLjA2MjkzOSAyMi43NSwxMzkgYyAyMi43NSwxNDQuOTM3MDYxIDE3LjkzNzA2MTEsMTQ5Ljc1IDEyLDE0OS43NSBjIDYuMDYyOTM4OTQsMTQ5Ljc1IDEuMjUsMTQ0LjkzNzA2MSAxLjI1LDEzOSB6XCI7XG52YXIgcG9pbnQ4ID0gXCJtIDk1LjI0MjkyMDksMTMzLjE0MTYyNCBjIDk1LjU1Mzc5MjcsMTI3LjIwOTgzNiA5MC43NTYyNjk4LDEyMi4xNDE3MzYgODQuNTMzNzQwNCwxMjEuODE1NjI3IGMgNzguMzExMjExLDEyMS40ODk1MTggNzMuMDEwMjA4NiwxMjYuMDI4Mzc3IDcyLjY5OTMzNjgsMTMxLjk2MDE2NSBjIDcyLjM4ODQ2NSwxMzcuODkxOTUzIDc3LjE4NTk4NzksMTQyLjk2MDA1MyA4My40MDg1MTczLDE0My4yODYxNjIgYyA4OS42MzEwNDY3LDE0My42MTIyNzEgOTQuOTMyMDQ5LDEzOS4wNzM0MTIgOTUuMjQyOTIwOSwxMzMuMTQxNjI0IHogbSA3My4xOTg2NTE2LDEzMS45ODYzMzMgYyA3My40OTQ3NzExLDEyNi4zMzYwMzYgNzguNTU1Mzg4LDEyMi4wMDMwMDEgODQuNTA3NTcyNCwxMjIuMzE0OTQyIGMgOTAuNDU5NzU2NywxMjIuNjI2ODgzIDk1LjAzOTcyNTYsMTI3LjQ2NTE1OSA5NC43NDM2MDYxLDEzMy4xMTU0NTYgYyA5NC40NDc0ODY2LDEzOC43NjU3NTQgODkuMzg2ODY5NiwxNDMuMDk4Nzg4IDgzLjQzNDY4NTMsMTQyLjc4Njg0NyBjIDc3LjQ4MjUwMDksMTQyLjQ3NDkwNyA3Mi45MDI1MzIsMTM3LjYzNjYzIDczLjE5ODY1MTYsMTMxLjk4NjMzMyB6XCI7XG52YXIgcG9pbnQ5ID0gXCJtIDE2Ny43MjgwMjMsMTM1LjA2MjQ3OSBjIDE2OC4wMzgzMjcsMTI5LjE0MTUzMyAxNjMuNDgyMTE3LDEyNC4wODk4NjQgMTU3LjU1MTY2MiwxMjMuNzc5MDYyIGMgMTUxLjYyMTIwOCwxMjMuNDY4MjYxIDE0Ni41NjE5MTQsMTI4LjAxNjAwMiAxNDYuMjUxNjEsMTMzLjkzNjk0OCBjIDE0NS45NDEzMDcsMTM5Ljg1Nzg5NCAxNTAuNDk3NTE3LDE0NC45MDk1NjIgMTU2LjQyNzk3MSwxNDUuMjIwMzY0IGMgMTYyLjM1ODQyNiwxNDUuNTMxMTY2IDE2Ny40MTc3MiwxNDAuOTgzNDI1IDE2Ny43MjgwMjMsMTM1LjA2MjQ3OSB6IG0gMTQ2Ljc1MDkyNSwxMzMuOTYzMTE2IGMgMTQ3LjA0Njc2NywxMjguMzE4MTIxIDE1MS44NzA2MTcsMTIzLjk4MjAxOCAxNTcuNTI1NDk0LDEyNC4yNzgzNzcgYyAxNjMuMTgwMzcyLDEyNC41NzQ3MzcgMTY3LjUyNDU1LDEyOS4zOTEzMTcgMTY3LjIyODcwOSwxMzUuMDM2MzExIGMgMTY2LjkzMjg2NywxNDAuNjgxMzA1IDE2Mi4xMDkwMTcsMTQ1LjAxNzQwOSAxNTYuNDU0MTM5LDE0NC43MjEwNSBjIDE1MC43OTkyNjIsMTQ0LjQyNDY5IDE0Ni40NTUwODMsMTM5LjYwODExIDE0Ni43NTA5MjUsMTMzLjk2MzExNiB6XCI7XG52YXIgbGV0dHJlc0FOVCA9IFwibSAyNS42MTg1NDE0LDMxLjk0NDcyNjYgbCAyNS42MTg1NDE0LDMxLjEwMDM5MDYgbCAyNS4zODgyNjgsMzAuNDA5NTcwMyBsIDI0LjkyNzcyMTEsMzAuMDI1NzgxMiBsIDI0LjU0MzkzMjEsMjkuNTY1MjM0NCBsIDI0LjA4MzM4NTIsMjkuMTgxNDQ1MyBsIDIzLjY5OTU5NjEsMjguNzIwODk4NCBsIDIzLjAwODc3NTgsMjguNDkwNjI1IGwgMjIuMTQ1MzY3OCwyOC40OTA2MjUgbCAyMS4xNTEyNzExLDI4LjQ5MDYyNSBsIDIwLjQ3NTc2OCwyOC40OTA2MjUgbCAxOS44NjE3MDU1LDI4LjcyMDg5ODQgbCAxOS40Nzc5MTY0LDI5LjE4MTQ0NTMgbCAxOC43ODcwOTYxLDI5LjMzNDk2MDkgbCAxNy41NDc4NTc5LDI5LjMzNDk2MDkgbCAxNi42Njk5MDcsMjkuMzM0OTYwOSBsIDE2LjAwMjQ5NzEsMjkuMzM0OTYwOSBsIDE1LjMzNDk4MDIsMjkuMzM0OTYwOSBsIDE0LjU2NTQxNjQsMjkuMzM0OTYwOSBsIDEzLjk1MTM1MzksMjkuNTY1MjM0NCBsIDEzLjU2NzU2NDksMzAuMDI1NzgxMiBsIDEyLjg3Njc0NDYsMzAuMTc5Mjk2OSBsIDEyLjAzMjQwODYsMzAuMTc5Mjk2OSBsIDExLjQxODM0NjEsMzAuNDA5NTcwMyBsIDExLjAzNDU1NzEsMzAuODcwMTE3MiBsIDEwLjM0MzczNjcsMzEuMDIzNjMyOCBsIDkuNzI5Njc0MjQsMzEuMjUzOTA2MiBsIDkuMzQ1ODg1MTgsMzEuNzE0NDUzMSBsIDguODg1MzM4MywzMi4wOTgyNDIyIGwgOC41MDE1NDkyNCwzMi41NTg3ODkxIGwgOC4wNDEwMDIzNywzMi45NDI1NzgxIGwgNy44ODc0ODY3NCwzMy42MzMzOTg0IGwgNy42NTcyMTMzLDM0LjI0NzQ2MDkgbCA3LjE5NjY2NjQzLDM0LjYzMTI1IGwgNy4wNDMxNTA4LDM1LjMyMjA3MDMgbCA3LjA0MzE1MDgsMzYuMTY2NDA2MiBsIDcuMTk2NjY2NDMsMzYuNzgwNDY4NyBsIDcuNjU3MjEzMywzNy4xNjQyNTc4IGwgOC4wNDEwMDIzNywzNy42MjQ4MDQ3IGwgOC41MDE1NDkyNCwzOC4wMDg1OTM3IGwgOC44ODUzMzgzLDM4LjQ2OTE0MDYgbCA5LjM0NTg4NTE4LDM4Ljg1MjkyOTcgbCA5LjcyOTY3NDI0LDM5LjMxMzQ3NjYgbCAxMC4zNDM3MzY3LDM5LjQ2Njk5MjIgbCAxMS40NzI1OTYxLDM5LjQ2Njk5MjIgbCAxMi4zMTIyOTYxLDM5LjQ2Njk5MjIgbCAxMi44NzQzODE4LDM5LjQ2Njk5MjIgbCAxMy42NDk2ODczLDM5LjQ2Njk5MjIgbCAxNC41Mzk4NTI4LDM5LjQ2Njk5MjIgbCAxNS40MDk3NTI0LDM5LjQ2Njk5MjIgbCAxNi4xMDA1NzI3LDM5LjMxMzQ3NjYgbCAxNi40ODQzNjE3LDM4Ljg1MjkyOTcgbCAxNy4wOTg0MjQyLDM4LjYyMjY1NjIgbCAxNy43ODkyNDQ2LDM4LjQ2OTE0MDYgbCAxOC4xNzMwMzM2LDM4LjAwODU5MzcgbCAxOC42MzM1ODA1LDM3LjYyNDgwNDcgbCAxOS4wMTczNjk2LDM3LjE2NDI1NzggbCAxOS42MzE0MzIxLDM2LjkzMzk4NDQgbCAyMC4zMjIyNTI0LDM2Ljc4MDQ2ODcgbCAyMC43MDYwNDE0LDM2LjMxOTkyMTkgbCAyMS4zMjAxMDM5LDM2LjA4OTY0ODQgbCAyMi4wMTA5MjQyLDM1LjkzNjEzMjggbCAyMi4zOTQ3MTMzLDM1LjQ3NTU4NTkgbCAyMy4wMDg3NzU4LDM1LjI0NTMxMjUgbCAyMy42OTk1OTYxLDM1LjA5MTc5NjkgbCAyMy45Mjk4Njk2LDM0LjQ3NzczNDQgbCAyNC4wODMzODUyLDMzLjc4NjkxNDEgbCAyNC41NDM5MzIxLDMzLjQwMzEyNSBsIDI0LjkyNzcyMTEsMzIuOTQyNTc4MSBsIDI1LjM4ODI2OCwzMi41NTg3ODkxIGwgMjUuNjE4NTQxNCwzMS45NDQ3MjY2IHogbSAzNy41OTI3NjAyLDQzLjkxODk0NTMgbCAzNy4yMDg5NzExLDQ0LjM3OTQ5MjIgbCAzNi43NDg0MjQyLDQ0Ljc2MzI4MTMgbCAzNi4zNjQ2MzUyLDQ1LjIyMzgyODEgbCAzNS42NzM4MTQ5LDQ1LjM3NzM0MzggbCAzNS4wNTk3NTI0LDQ1LjYwNzYxNzIgbCAzNC42NzU5NjMzLDQ2LjA2ODE2NDEgbCAzMy45ODUxNDMsNDYuMjIxNjc5NyBsIDMzLjM5NTE3Miw0Ni4yMjE2Nzk3IGwgMzIuNDg3NDIxOSw0Ni4yMjE2Nzk3IGwgMzEuNDUyMTM1Miw0Ni4yMjE2Nzk3IGwgMzAuODM4MDcyNyw0Ni4wNjgxNjQxIGwgMzAuNDU0MjgzNiw0NS42MDc2MTcyIGwgMjkuNzYzNDYzMyw0NS4zNzczNDM4IGwgMjkuMTQ5NDAwOCw0NS4yMjM4MjgxIGwgMjguNzY1NjExNyw0NC43NjMyODEzIGwgMjguMzA1MDY0OSw0NC4zNzk0OTIyIGwgMjcuOTIxMjc1OCw0My45MTg5NDUzIGwgMjcuMjMwNDU1NSw0My42ODg2NzE5IGwgMjYuNjE2MzkzLDQzLjUzNTE1NjMgbCAyNi4yMzI2MDM5LDQzLjA3NDYwOTQgbCAyNS43NzIwNTcxLDQzLjA3NDYwOTQgbCAyNS4zODgyNjgsNDMuNTM1MTU2MyBsIDI0LjY5NzQ0NzcsNDMuNjg4NjcxOSBsIDIzLjg1MzExMTcsNDMuNjg4NjcxOSBsIDIzLjIzOTA0OTIsNDMuOTE4OTQ1MyBsIDIyLjg1NTI2MDIsNDQuMzc5NDkyMiBsIDIyLjE2NDQzOTksNDQuNTMzMDA3OCBsIDIxLjU1MDM3NzQsNDQuNzYzMjgxMyBsIDIxLjE2NjU4ODMsNDUuMjIzODI4MSBsIDIwLjQ3NTc2OCw0NS4zNzczNDM4IGwgMTkuNjMxNDMyMSw0NS4zNzczNDM4IGwgMTkuMDE3MzY5Niw0NS42MDc2MTcyIGwgMTguNjMzNTgwNSw0Ni4wNjgxNjQxIGwgMTcuOTQyNzYwMiw0Ni4yMjE2Nzk3IGwgMTcuMTI3OTA4Nyw0Ni4yMjE2Nzk3IGwgMTYuMjU0MDg4Myw0Ni4yMjE2Nzk3IGwgMTUuNjQwMDI1OCw0Ni40NTE5NTMxIGwgMTUuMjU2MjM2Nyw0Ni45MTI1IGwgMTQuNzk1Njg5OSw0Ni45MTI1IGwgMTQuNDExOTAwOCw0Ni40NTE5NTMxIGwgMTMuNzIxMDgwNSw0Ni4yMjE2Nzk3IGwgMTIuNzkyNjgzNyw0Ni4yMjE2Nzk3IGwgMTIuMDE2OTgsNDYuMjIxNjc5NyBsIDExLjM2NDY5LDQ2LjIyMTY3OTcgbCAxMC40ODQ5MSw0Ni4yMjE2Nzk3IGwgOS40OTk0MDA4LDQ2LjIyMTY3OTcgbCA4Ljg4NTMzODMsNDYuMDY4MTY0MSBsIDguNTAxNTQ5MjQsNDUuNjA3NjE3MiBsIDcuODEwNzI4OTMsNDUuMzc3MzQzOCBsIDcuMTk2NjY2NDMsNDUuMjIzODI4MSBsIDYuODEyODc3MzcsNDQuNzYzMjgxMyBsIDYuMTIyMDU3MDUsNDQuNTMzMDA3OCBsIDUuNTA3OTk0NTUsNDQuMzc5NDkyMiBsIDUuMTI0MjA1NDksNDMuOTE4OTQ1MyBsIDQuNjYzNjU4NjIsNDMuNTM1MTU2MyBsIDQuMjc5ODY5NTUsNDMuMDc0NjA5NCBsIDMuODE5MzIyNjgsNDIuNjkwODIwMyBsIDMuNDM1NTMzNjIsNDIuMjMwMjczNCBsIDIuOTc0OTg2NzQsNDEuODQ2NDg0NCBsIDIuNTkxMTk3NjgsNDEuMzg1OTM3NSBsIDIuMTMwNjUwOCw0MS4wMDIxNDg0IGwgMS45NzcxMzUxOCw0MC4zODgwODU5IGwgMS43NDY4NjE3NCwzOS42OTcyNjU2IGwgMS4yODYzMTQ4NywzOS4zMTM0NzY2IGwgMS4xMzI3OTkyNCwzOC42OTk0MTQxIGwgMS4xMzI3OTkyNCwzNy44NTUwNzgxIGwgMC45MDI1MjU4MDQsMzcuMTY0MjU3OCBsIDAuNDQxOTc4OTI5LDM2Ljc4MDQ2ODggbCAwLjI4ODQ2MzMwNCwzNi4xNjY0MDYzIGwgMC4yODg0NjMzMDQsMzUuMjUyMTQzMyBsIDAuMjg4NDYzMzA0LDM0LjM5NzA4NzUgbCAwLjI4ODQ2MzMwNCwzMy42MzMzOTg0IGMgMC4yODg0NjMzMDQsMzMuNTAwODcwNSAwLjQ0MTk3ODkyOSwzMi45NDI1NzgxIDAuNDQxOTc4OTI5LDMyLjk0MjU3ODEgbCAwLjkwMjUyNTgwNCwzMi41NTg3ODkxIGwgMS4xMzI3OTkyNCwzMS45NDQ3MjY2IGwgMS4xMzI3OTkyNCwzMS4xMDAzOTA2IGwgMS4yODYzMTQ4NywzMC40MDk1NzAzIGwgMS43NDY4NjE3NCwzMC4wMjU3ODEzIGwgMS45NzcxMzUxOCwyOS40MTE3MTg4IGwgMi4xMzA2NTA4LDI4LjcyMDg5ODQgbCAyLjU5MTE5NzY4LDI4LjMzNzEwOTQgbCAyLjk3NDk4Njc0LDI3Ljg3NjU2MjUgbCAzLjQzNTUzMzYyLDI3LjQ5Mjc3MzQgbCAzLjgxOTMyMjY4LDI3LjAzMjIyNjYgbCA0LjI3OTg2OTU1LDI2LjY0ODQzNzUgbCA0LjY2MzY1ODYyLDI2LjE4Nzg5MDYgbCA1LjI3NzcyMTEyLDI1Ljk1NzYxNzIgbCA2LjEyMjA1NzA1LDI1Ljk1NzYxNzIgbCA2LjgxMjg3NzM3LDI1LjgwNDEwMTYgbCA3LjE5NjY2NjQzLDI1LjM0MzU1NDcgbCA3LjgxMDcyODkzLDI1LjExMzI4MTMgbCA4LjY1NTA2NDg3LDI1LjExMzI4MTMgbCA5LjM0NTg4NTE4LDI0Ljk1OTc2NTYgbCA5LjcyOTY3NDI0LDI0LjQ5OTIxODggbCAxMC4zNDM3MzY3LDI0LjI2ODk0NTMgbCAxMC45MzcwNTg1LDI0LjI2ODk0NTMgbCAxMS42Njk1Mjk0LDI0LjI2ODk0NTMgbCAxMi4xMjQ0NDA0LDI0LjI2ODk0NTMgbCAxMi43ODkwOTI1LDI0LjI2ODk0NTMgbCAxMy43MjEwODA1LDI0LjI2ODk0NTMgYyAxMy45NTY2MDg0LDI0LjI2ODk0NTMgMTQuNDExOTAwOCwyNC4xMTU0Mjk3IDE0LjQxMTkwMDgsMjQuMTE1NDI5NyBsIDE0Ljc5NTY4OTksMjMuNjU0ODgyOCBsIDE1LjQwOTc1MjQsMjMuNDI0NjA5NCBsIDE1Ljg5MTQ5MjksMjMuNDI0NjA5NCBsIDE2LjM3ODAzNDIsMjMuNDI0NjA5NCBsIDE2LjcxNjUwNjgsMjMuNDI0NjA5NCBsIDE3LjMyNjQ2MTgsMjMuNDI0NjA5NCBsIDE3Ljk0Mjc2MDIsMjMuNDI0NjA5NCBjIDE4LjE3ODI4ODEsMjMuNDI0NjA5NCAxOC42MzM1ODA1LDIzLjI3MTA5MzggMTguNjMzNTgwNSwyMy4yNzEwOTM4IGwgMTkuMDE3MzY5NiwyMi44MTA1NDY5IGwgMTkuNjMxNDMyMSwyMi41ODAyNzM0IGwgMjAuMTI1ODg3NSwyMi41ODAyNzM0IGwgMjAuNjU1NjYxNiwyMi41ODAyNzM0IGwgMjEuMzAyNjMzOCwyMi41ODAyNzM0IGwgMjIuNDkxNjczMiwyMi41ODAyNzM0IGwgMjMuNDM0NzgzMSwyMi41ODAyNzM0IGwgMjMuODUzMTExNywyMi41ODAyNzM0IGwgMjQuNTQzOTMyMSwyMi40MjY3NTc4IGwgMjQuOTI3NzIxMSwyMS45NjYyMTA5IGwgMjUuMzg4MjY4LDIxLjU4MjQyMTkgbCAyNS43NzIwNTcxLDIxLjEyMTg3NSBsIDI2LjM4NjExOTYsMjAuODkxNjAxNiBsIDI3LjA3NjkzOTksMjAuNzM4MDg1OSBsIDI3LjMwNzIxMzMsMjAuMTI0MDIzNCBsIDI3LjMwNzIxMzMsMTkuMzIwMTA2NyBsIDI3LjMwNzIxMzMsMTguNDM1MzUxNiBsIDI3LjA3NjkzOTksMTcuNzQ0NTMxMyBsIDI2LjYxNjM5MywxNy4zNjA3NDIyIGwgMjYuNjE2MzkzLDE2LjkwMDE5NTMgbCAyNy4wNzY5Mzk5LDE2LjUxNjQwNjMgbCAyNy4wNzY5Mzk5LDE2LjA1NTg1OTQgbCAyNi42MTYzOTMsMTUuNjcyMDcwMyBsIDI2LjQ2Mjg3NzQsMTUuMDU4MDA3OCBsIDI2LjIzMjYwMzksMTQuMzY3MTg3NSBsIDI1Ljc3MjA1NzEsMTMuOTgzMzk4NCBsIDI1LjM4ODI2OCwxMy41MjI4NTE2IGwgMjQuOTI3NzIxMSwxMy4xMzkwNjI1IGwgMjQuNTQzOTMyMSwxMi42Nzg1MTU2IGwgMjMuODUzMTExNywxMi40NDgyNDIyIGwgMjMuMDA4Nzc1OCwxMi40NDgyNDIyIGwgMjIuMzk0NzEzMywxMi4yOTQ3MjY2IGwgMjIuMDEwOTI0MiwxMS44MzQxNzk3IGwgMjEuNTUwMzc3NCwxMS44MzQxNzk3IGwgMjEuMTY2NTg4MywxMi4yOTQ3MjY2IGwgMjAuNDc1NzY4LDEyLjQ0ODI0MjIgbCAxOS41MzU2MTUxLDEyLjQ0ODI0MjIgbCAxOC42Nzg5MDE4LDEyLjQ0ODI0MjIgbCAxOC4xMDA2NDYxLDEyLjQ0ODI0MjIgbCAxNy40OTA2OTExLDEyLjQ0ODI0MjIgbCAxNi44NjkyMDI3LDEyLjQ0ODI0MjIgbCAxNi4yOTYyNjQ5LDEyLjQ0ODI0MjIgbCAxNS43NTU4NTUyLDEyLjQ0ODI0MjIgbCAxNS4xNjQ0MDg4LDEyLjQ0ODI0MjIgbCAxNC41NjU0MTY0LDEyLjQ0ODI0MjIgbCAxMy45NTEzNTM5LDEyLjY3ODUxNTYgbCAxMy41Njc1NjQ5LDEzLjEzOTA2MjUgbCAxMy4xMDcwMTgsMTMuNTIyODUxNiBsIDEyLjcyMzIyODksMTMuOTgzMzk4NCBsIDEyLjI2MjY4MjEsMTQuMzY3MTg3NSBsIDEyLjEwOTE2NjQsMTUuMDU4MDA3OCBsIDEyLjEwOTE2NjQsMTUuNzAyMDg5OSBsIDEyLjEwOTE2NjQsMTYuNjk0NTc5IGwgMTIuMTA5MTY2NCwxNy41OTEwMTU2IGwgMTEuODc4ODkzLDE4LjIwNTA3ODEgbCAxMS40MTgzNDYxLDE4LjU4ODg2NzIgbCAxMS4wMzQ1NTcxLDE5LjA0OTQxNDEgbCAxMC41NzQwMTAyLDE5LjQzMzIwMzEgbCAxMC4xOTAyMjExLDE5Ljg5Mzc1IGwgOS40OTk0MDA4LDIwLjA0NzI2NTYgbCA4LjY0NDQ3ODQsMjAuMDQ3MjY1NiBsIDcuODEwNzI4OTMsMjAuMDQ3MjY1NiBsIDcuMTk2NjY2NDMsMTkuODkzNzUgbCA2LjgxMjg3NzM3LDE5LjQzMzIwMzEgbCA2LjEyMjA1NzA1LDE5LjIwMjkyOTcgbCA1LjUwNzk5NDU1LDE5LjA0OTQxNDEgbCA1LjEyNDIwNTQ5LDE4LjU4ODg2NzIgbCA0LjY2MzY1ODYyLDE4LjIwNTA3ODEgbCA0LjUxMDE0Mjk5LDE3LjU5MTAxNTYgbCA0LjUxMDE0Mjk5LDE2Ljc0NjY3OTcgbCA0LjI3OTg2OTU1LDE2LjA1NTg1OTQgbCAzLjgxOTMyMjY4LDE1LjY3MjA3MDMgbCAzLjY2NTgwNzA1LDE1LjA1ODAwNzggbCAzLjY2NTgwNzA1LDE0LjIxMzY3MTkgbCAzLjgxOTMyMjY4LDEzLjUyMjg1MTYgbCA0LjI3OTg2OTU1LDEzLjEzOTA2MjUgbCA0LjUxMDE0Mjk5LDEyLjUyNSBsIDQuNjYzNjU4NjIsMTEuODM0MTc5NyBsIDUuMTI0MjA1NDksMTEuNDUwMzkwNiBsIDUuNTA3OTk0NTUsMTAuOTg5ODQzOCBsIDUuOTY4NTQxNDMsMTAuNjA2MDU0NyBsIDYuMzUyMzMwNDksMTAuMTQ1NTA3OCBsIDYuOTY2MzkyOTksOS45MTUyMzQzOCBsIDcuNjU3MjEzMyw5Ljc2MTcxODc1IGwgOC4wNDEwMDIzNyw5LjMwMTE3MTg4IGwgOC42NTUwNjQ4Nyw5LjA3MDg5ODQ0IGwgOS4zNDU4ODUxOCw4LjkxNzM4MjgxIGwgOS43Mjk2NzQyNCw4LjQ1NjgzNTk0IGwgMTAuMzQzNzM2Nyw4LjIyNjU2MjUgbCAxMC45ODkwNjIxLDguMjI2NTYyNSBsIDExLjQzNjEwMDEsOC4yMjY1NjI1IGwgMTEuOTQyODA3NSw4LjIyNjU2MjUgbCAxMi4yNjgxNTgzLDguMjI2NTYyNSBsIDEyLjk1NDk3OTMsOC4yMjY1NjI1IGwgMTMuNzIxMDgwNSw4LjIyNjU2MjUgbCAxNC40MTE5MDA4LDguMDczMDQ2ODggbCAxNC43OTU2ODk5LDcuNjEyNSBsIDE1LjI1NjIzNjcsNy42MTI1IGwgMTUuNjQwMDI1OCw4LjA3MzA0Njg4IGwgMTYuMjU0MDg4Myw4LjIyNjU2MjUgbCAxNi45NjA5MTcsOC4yMjY1NjI1IGwgMTcuNDkxNTE5OSw4LjIyNjU2MjUgbCAxOC4wMjkwOTgsOC4yMjY1NjI1IGwgMTguNzkzMjY4Myw4LjIyNjU2MjUgbCAxOS42MzE0MzIxLDguMjI2NTYyNSBsIDIwLjMyMjI1MjQsOC40NTY4MzU5NCBsIDIwLjcwNjA0MTQsOC45MTczODI4MSBsIDIxLjMyMDEwMzksOS4wNzA4OTg0NCBsIDIyLjA0NDA4MjgsOS4wNzA4OTg0NCBsIDIyLjg3OTQ1Niw5LjA3MDg5ODQ0IGwgMjMuODUzMTExNyw5LjA3MDg5ODQ0IGMgMjQuMDk1MDQ3NCw5LjA3MDg5ODQ0IDI0LjU0MzkzMjEsOS4zMDExNzE4OCAyNC41NDM5MzIxLDkuMzAxMTcxODggbCAyNC45Mjc3MjExLDkuNzYxNzE4NzUgbCAyNS41NDE3ODM2LDkuOTE1MjM0MzggbCAyNi4yMzI2MDM5LDEwLjE0NTUwNzggbCAyNi42MTYzOTMsMTAuNjA2MDU0NyBsIDI3LjA3NjkzOTksMTAuOTg5ODQzOCBsIDI3LjQ2MDcyODksMTEuNDUwMzkwNiBsIDI3LjkyMTI3NTgsMTEuODM0MTc5NyBsIDI4LjMwNTA2NDksMTIuMjk0NzI2NiBsIDI4Ljc2NTYxMTcsMTIuNjc4NTE1NiBsIDI5LjE0OTQwMDgsMTMuMTM5MDYyNSBsIDI5LjYwOTk0NzcsMTMuNTIyODUxNiBsIDI5Ljg0MDIyMTEsMTQuMjEzNjcxOSBsIDI5Ljk5MzczNjcsMTQuODI3NzM0NCBsIDMwLjQ1NDI4MzYsMTUuMjExNTIzNCBsIDMwLjY4NDU1NzEsMTUuOTAyMzQzOCBsIDMwLjY4NDU1NzEsMTYuNTAzMjc3NCBjIDMwLjY4NDU1NzEsMTYuNzI3MTUzMiAzMC42ODQ1NTcxLDE2Ljk1MTAyOSAzMC42ODQ1NTcxLDE3LjE3NDkwNDggbCAzMC42ODQ1NTcxLDE3Ljg3NjU3NDEgbCAzMC42ODQ1NTcxLDE4LjU0MjA1NDkgbCAzMC42ODQ1NTcxLDE5LjE2NTIwMDggYyAzMC42ODQ1NTcxLDE5LjcwOTExNDcgMzAuNjg0NTU3MSwyMC4xMjQwMjM0IDMwLjY4NDU1NzEsMjAuMTI0MDIzNCBsIDMwLjgzODA3MjcsMjAuNzM4MDg1OSBsIDMxLjI5ODYxOTYsMjEuMTIxODc1IGwgMzEuNTI4ODkzLDIxLjgxMjY5NTMgbCAzMS41Mjg4OTMsMjIuNjA2MDczNyBjIDMxLjUyODg5MywyMi44NjI4NjkgMzEuNTI4ODkzLDIzLjExOTY2NDMgMzEuNTI4ODkzLDIzLjM3NjQ1OTYgbCAzMS41Mjg4OTMsMjQuMTEwNTg4IGMgMzEuNTI4ODkzLDI0LjM3NzA3NDIgMzEuNTI4ODkzLDI0LjY1NTU4NDYgMzEuNTI4ODkzLDI0Ljk0MzQ3NSBsIDMxLjUyODg5MywyNS44ODc0MTM2IGwgMzEuNTI4ODkzLDI2Ljc3NDE2ODkgYyAzMS41Mjg4OTMsMjcuMDc0NDA0MSAzMS41Mjg4OTMsMjcuMzc0NjM5NCAzMS41Mjg4OTMsMjcuNjc0ODc0NyBsIDMxLjUyODg5MywyOC41ODcxMTM4IGwgMzEuNTI4ODkzLDI5LjI2MjEyNTIgbCAzMS41Mjg4OTMsMzAuMDE5MjUxMiBjIDMxLjUyODg5MywzMC4yNjQ0NzY4IDMxLjUyODg5MywzMC41MTI4MDYzIDMxLjUyODg5MywzMC43NjEyNTI3IGwgMzEuNTI4ODkzLDMxLjU5MTY1MzQgbCAzMS41Mjg4OTMsMzIuMjc1NjQyOSBsIDMxLjUyODg5MywzMy4xNTAwMzYxIGwgMzEuNTI4ODkzLDMzLjgxMzg1OTQgbCAzMS41Mjg4OTMsMzQuNTUzMzc0NiBjIDMxLjUyODg5MywzNS4wMzEwNzc4IDMxLjUyODg5MywzNS4zMjIwNzAzIDMxLjUyODg5MywzNS4zMjIwNzAzIGwgMzEuNjgyNDA4NiwzNS45MzYxMzI4IGwgMzIuMTQyOTU1NSwzNi4zMTk5MjE5IGwgMzIuNTI2NzQ0NiwzNi43ODA0Njg4IGwgMzIuOTg3MjkxNCwzNy4xNjQyNTc4IGwgMzMuMzcxMDgwNSwzNy42MjQ4MDQ3IGwgMzMuOTg1MTQzLDM3Ljc3ODMyMDMgbCAzNC42NzU5NjMzLDM4LjAwODU5MzggbCAzNS4wNTk3NTI0LDM4LjQ2OTE0MDYgbCAzNS42NzM4MTQ5LDM4LjYyMjY1NjMgbCAzNi4zNjQ2MzUyLDM4Ljg1MjkyOTcgbCAzNi43NDg0MjQyLDM5LjMxMzQ3NjYgbCAzNy4zNjI0ODY3LDM5LjQ2Njk5MjIgbCAzOC4wNTMzMDcxLDM5LjY5NzI2NTYgbCAzOC4yODM1ODA1LDQwLjM4ODA4NTkgbCAzOC4yODM1ODA1LDQxLjAyMTIzMjMgbCAzOC4yODM1ODA1LDQyLjA3NTM5MzcgbCAzOC4yODM1ODA1LDQyLjkyMTA5MzggbCAzOC4wNTMzMDcxLDQzLjUzNTE1NjMgbCAzNy41OTI3NjAyLDQzLjkxODk0NTMgeiBtIDg3LjQxODE4NTcsNDAuNTQxNjAxNiBsIDg3LjI2NDY3MDEsNDEuMjMyNDIxOSBsIDg3LjAzNDM5NjYsNDEuODQ2NDg0NCBsIDg2LjU3Mzg0OTcsNDIuMjMwMjczNCBsIDg2LjE5MDA2MDcsNDIuNjkwODIwMyBsIDg1LjcyOTUxMzgsNDMuMDc0NjA5NCBsIDg1LjM0NTcyNDcsNDMuNTM1MTU2MyBsIDg0LjY1NDkwNDQsNDMuNjg4NjcxOSBsIDg0LjA0MDg0MTksNDMuOTE4OTQ1MyBsIDgzLjY1NzA1MjksNDQuMzc5NDkyMiBsIDgyLjk2NjIzMjYsNDQuNTMzMDA3OCBsIDgyLjM4MjgxOTEsNDQuNTMzMDA3OCBsIDgxLjc2MTU0NjMsNDQuNTMzMDA3OCBsIDgxLjE0OTM5OTMsNDQuNTMzMDA3OCBsIDgwLjM5OTExMjIsNDQuNTMzMDA3OCBsIDc5Ljc2OTg3NjcsNDQuNTMzMDA3OCBsIDc5LjI2NzQxODgsNDQuNTMzMDA3OCBsIDc4LjU5OTI2NDIsNDQuNTMzMDA3OCBsIDc3Ljk0NzEyNDYsNDQuNTMzMDA3OCBsIDc3LjI1OTczNDIsNDQuNTMzMDA3OCBsIDc2LjY2OTQxNzcsNDQuNTMzMDA3OCBsIDc2LjE1ODk5Nyw0NC41MzMwMDc4IGwgNzUuNTg0NjA2LDQ0LjUzMzAwNzggbCA3NS4wMzQxOTI2LDQ0LjUzMzAwNzggbCA3NC41MDc3NTcsNDQuNTMzMDA3OCBsIDc0LjI3MTIwMDksNDQuNTMzMDA3OCBsIDczLjcwOTYwMzksNDQuNTMzMDA3OCBsIDcyLjgzNDIwMTMsNDQuNTMzMDA3OCBjIDcyLjg2ODY0MSw0NC41MDAwNDEzIDcyLjIyMDEzODgsNDQuMzc5NDkyMiA3Mi4yMjAxMzg4LDQ0LjM3OTQ5MjIgbCA3MS44MzYzNDk3LDQzLjkxODk0NTMgbCA3MS4zNzU4MDI5LDQzLjUzNTE1NjMgbCA3MS4yMjIyODcyLDQyLjkyMTA5MzggbCA3MC45OTIwMTM4LDQyLjIzMDI3MzQgbCA3MC41MzE0NjY5LDQxLjg0NjQ4NDQgbCA3MC4zNzc5NTEzLDQxLjIzMjQyMTkgbCA3MC41MzE0NjY5LDQwLjU0MTYwMTYgbCA3MC45OTIwMTM4LDQwLjE1NzgxMjUgbCA3MS4zNzU4MDI5LDM5LjY5NzI2NTYgbCA3MS44MzYzNDk3LDM5LjMxMzQ3NjYgbCA3Mi4yMjAxMzg4LDM4Ljg1MjkyOTcgbCA3Mi42ODA2ODU3LDM4LjQ2OTE0MDYgbCA3My4wNjQ0NzQ3LDM4LjAwODU5MzggbCA3My42Nzg1MzcyLDM3Ljc3ODMyMDMgbCA3NC4zNjkzNTc2LDM3LjYyNDgwNDcgbCA3NC43NTMxNDY2LDM3LjE2NDI1NzggbCA3NS4yMTM2OTM1LDM2Ljc4MDQ2ODggbCA3NS41OTc0ODI2LDM2LjMxOTkyMTkgbCA3Ni4wNTgwMjk0LDM1LjkzNjEzMjggbCA3Ni4yODgzMDI5LDM1LjMyMjA3MDMgbCA3Ni4yODgzMDI5LDM0LjQ3NzczNDQgbCA3Ni40NDE4MTg1LDMzLjc4NjkxNDEgbCA3Ni45MDIzNjU0LDMzLjQwMzEyNSBsIDc3LjEzMjYzODgsMzIuNzg5MDYyNSBsIDc3LjEzMjYzODgsMzIuMDU4NjUxNyBsIDc3LjEzMjYzODgsMzEuMTM0Nzk0NyBsIDc3LjEzMjYzODgsMzAuMDY4MDU1NyBsIDc3LjEzMjYzODgsMjkuMjYzMDEzNiBsIDc3LjEzMjYzODgsMjguMzg3MDIyNSBsIDc3LjEzMjYzODgsMjcuNTA0OTQ3NSBsIDc3LjEzMjYzODgsMjYuNjAyODMxNSBsIDc3LjEzMjYzODgsMjUuNzcxNjY0NCBsIDc3LjEzMjYzODgsMjQuODIzNzQwMyBsIDc3LjEzMjYzODgsMjMuOTIxNzEzOCBsIDc3LjEzMjYzODgsMjMuMDkwNDU3MiBsIDc3LjEzMjYzODgsMjIuMTMxMzQ5NCBsIDc3LjEzMjYzODgsMjEuMTY0Mjc4OSBsIDc3LjEzMjYzODgsMjAuMjM4Mjc0NyBsIDc3LjEzMjYzODgsMTkuNTQyOTIxNSBsIDc3LjEzMjYzODgsMTguOTIxNjQ4OCBsIDc3LjEzMjYzODgsMTguMTc5NDEzOSBsIDc3LjEzMjYzODgsMTcuNTkxMDE1NiBsIDc3LjI4NjE1NDQsMTYuOTAwMTk1MyBsIDc3Ljc0NjcwMTMsMTYuNTE2NDA2MyBsIDc3Ljk3Njk3NDcsMTUuOTAyMzQzOCBsIDc3Ljc0NjcwMTMsMTUuMjExNTIzNCBsIDc3LjI4NjE1NDQsMTQuODI3NzM0NCBsIDc3LjEzMjYzODgsMTQuMjEzNjcxOSBsIDc2LjkwMjM2NTQsMTMuNTIyODUxNiBsIDc2LjQ0MTgxODUsMTMuMTM5MDYyNSBsIDc2LjA1ODAyOTQsMTIuNjc4NTE1NiBsIDc1LjM2NzIwOTEsMTIuNDQ4MjQyMiBsIDc0Ljc1MzE0NjYsMTIuMjk0NzI2NiBsIDc0LjM2OTM1NzYsMTEuODM0MTc5NyBsIDczLjY3ODUzNzIsMTEuNjAzOTA2MyBsIDczLjA2NDQ3NDcsMTEuNDUwMzkwNiBsIDcyLjY4MDY4NTcsMTAuOTg5ODQzOCBsIDcxLjk4OTg2NTQsMTAuNzU5NTcwMyBsIDcxLjIxMzUwODEsMTAuNzU5NTcwMyBsIDcwLjI4ODQ4OCwxMC43NTk1NzAzIGwgNjkuNzcxMTc4MiwxMC43NTk1NzAzIGwgNjkuMDUzOTk0NiwxMC43NTk1NzAzIGwgNjguMjk2ODE4NCwxMC43NTk1NzAzIGwgNjcuNTY1ODU2NiwxMC43NTk1NzAzIGwgNjYuOTA0NTkxMiwxMC43NTk1NzAzIGwgNjYuMDc5NTEzOCwxMC43NTk1NzAzIGMgNjYuMTE2Mzg5MywxMC44MTAyMDEgNjUuNDY1NDUxMywxMC45ODk4NDM4IDY1LjQ2NTQ1MTMsMTAuOTg5ODQzOCBsIDY1LjA4MTY2MjIsMTEuNDUwMzkwNiBsIDY0LjM5MDg0MTksMTEuNjAzOTA2MyBsIDYzLjc3Njc3OTQsMTEuODM0MTc5NyBsIDYzLjYyMzI2MzgsMTIuNTI1IGwgNjMuMzkyOTkwNCwxMy4xMzkwNjI1IGwgNjIuOTMyNDQzNSwxMy41MjI4NTE2IGwgNjIuNTQ4NjU0NCwxMy45ODMzOTg0IGwgNjIuMDg4MTA3NiwxNC4zNjcxODc1IGwgNjEuNzA0MzE4NSwxNC44Mjc3MzQ0IGwgNjEuMjQzNzcxNiwxNS4yMTE1MjM0IGwgNjEuMDkwMjU2LDE1LjkwMjM0MzggbCA2MC44NTk5ODI2LDE2LjUxNjQwNjMgbCA2MC4zOTk0MzU3LDE2LjkwMDE5NTMgbCA2MC4wMTU2NDY2LDE3LjM2MDc0MjIgbCA1OS41NTUwOTk3LDE3Ljc0NDUzMTMgbCA1OS40MDE1ODQxLDE4LjQzNTM1MTYgbCA1OS4xNzEzMTA3LDE5LjA0OTQxNDEgbCA1OC43MTA3NjM4LDE5LjQzMzIwMzEgbCA1OC41NTcyNDgyLDIwLjEyNDAyMzQgbCA1OC4zMjY5NzQ3LDIwLjczODA4NTkgbCA1Ny44NjY0Mjc5LDIxLjEyMTg3NSBsIDU3LjcxMjkxMjIsMjEuODEyNjk1MyBsIDU3LjQ4MjYzODgsMjIuNDI2NzU3OCBsIDU3LjAyMjA5MTksMjIuODEwNTQ2OSBsIDU2Ljg2ODU3NjMsMjMuNTAxMzY3MiBsIDU2Ljg2ODU3NjMsMjMuOTg5NzEwMiBsIDU2Ljg2ODU3NjMsMjQuNjk4OTMxIGwgNTYuODY4NTc2MywyNS42NTQ4MTc5IGwgNTYuODY4NTc2MywyNi4zMjUxMTk3IGwgNTYuODY4NTc2MywyNy4wOTAyNTg3IGwgNTYuODY4NTc2MywyNy44MDA0NjM3IGwgNTYuODY4NTc2MywyOC41NjczODI4IGwgNTcuMDIyMDkxOSwyOS4xODE0NDUzIGwgNTcuNDgyNjM4OCwyOS41NjUyMzQ0IGwgNTcuNzEyOTEyMiwzMC4yNTYwNTQ3IGwgNTcuNzEyOTEyMiwzMS4xMDAzOTA2IGwgNTcuNDgyNjM4OCwzMS43MTQ0NTMxIGwgNTcuMDIyMDkxOSwzMi4wOTgyNDIyIGwgNTYuODY4NTc2MywzMi43ODkwNjI1IGwgNTcuMDIyMDkxOSwzMy40MDMxMjUgbCA1Ny40ODI2Mzg4LDMzLjc4NjkxNDEgbCA1Ny43MTI5MTIyLDM0LjQ3NzczNDQgbCA1Ny43MTI5MTIyLDM1LjMyMjA3MDMgbCA1Ny44NjY0Mjc5LDM1LjkzNjEzMjggbCA1OC4zMjY5NzQ3LDM2LjMxOTkyMTkgbCA1OC43MTA3NjM4LDM2Ljc4MDQ2ODggbCA1OS4xNzEzMTA3LDM3LjE2NDI1NzggbCA1OS41NTUwOTk3LDM3LjYyNDgwNDcgbCA2MC4xNjkxNjIyLDM3Ljc3ODMyMDMgbCA2MC44NTk5ODI2LDM4LjAwODU5MzggbCA2MS4yNDM3NzE2LDM4LjQ2OTE0MDYgbCA2MS43MDQzMTg1LDM4Ljg1MjkyOTcgbCA2Mi4wODgxMDc2LDM5LjMxMzQ3NjYgbCA2Mi41NDg2NTQ0LDM5LjY5NzI2NTYgbCA2Mi43Nzg5Mjc5LDQwLjM4ODA4NTkgbCA2Mi45MzI0NDM1LDQxLjAwMjE0ODQgbCA2My4zOTI5OTA0LDQxLjM4NTkzNzUgbCA2My4zOTI5OTA0LDQxLjg0NjQ4NDQgbCA2Mi45MzI0NDM1LDQyLjIzMDI3MzQgbCA2Mi43Nzg5Mjc5LDQyLjkyMTA5MzggbCA2Mi41NDg2NTQ0LDQzLjUzNTE1NjMgbCA2MS44NTc4MzQxLDQzLjY4ODY3MTkgbCA2MS4yNDM3NzE2LDQzLjkxODk0NTMgbCA2MC44NTk5ODI2LDQ0LjM3OTQ5MjIgbCA2MC4xNjkxNjIyLDQ0LjUzMzAwNzggbCA1OS4zMjQ4MjYzLDQ0LjUzMzAwNzggbCA1OC43MTA3NjM4LDQ0Ljc2MzI4MTMgbCA1OC4zMjY5NzQ3LDQ1LjIyMzgyODEgbCA1Ny44NjY0Mjc5LDQ1LjIyMzgyODEgbCA1Ny40ODI2Mzg4LDQ0Ljc2MzI4MTMgbCA1Ny4wMjIwOTE5LDQ0Ljc2MzI4MTMgbCA1Ni42MzgzMDI5LDQ1LjIyMzgyODEgbCA1NS45NDc0ODI2LDQ1LjM3NzM0MzggbCA1NS4xMDMxNDY2LDQ1LjM3NzM0MzggbCA1NC40ODkwODQxLDQ1LjIyMzgyODEgbCA1NC4xMDUyOTUxLDQ0Ljc2MzI4MTMgbCA1My42NDQ3NDgyLDQ0Ljc2MzI4MTMgbCA1My4yNjA5NTkxLDQ1LjIyMzgyODEgbCA1Mi41NzAxMzg4LDQ1LjM3NzM0MzggbCA1MS40OTc0OTExLDQ1LjM3NzM0MzggbCA1MC41OTA2MzMyLDQ1LjM3NzM0MzggbCA0OS45ODUzNzU0LDQ1LjM3NzM0MzggbCA0OS4xOTI3OTUxLDQ1LjM3NzM0MzggYyA0OS4yMjY4NDAyLDQ1LjM0NDMwMzkgNDguNTc4NzMyNiw0NS4yMjM4MjgxIDQ4LjU3ODczMjYsNDUuMjIzODI4MSBsIDQ4LjE5NDk0MzUsNDQuNzYzMjgxMyBsIDQ3LjUwNDEyMzIsNDQuNTMzMDA3OCBsIDQ2Ljg5MDA2MDcsNDQuMzc5NDkyMiBsIDQ2LjUwNjI3MTYsNDMuOTE4OTQ1MyBsIDQ2LjA0NTcyNDcsNDMuNTM1MTU2MyBsIDQ1LjY2MTkzNTcsNDMuMDc0NjA5NCBsIDQ1LjIwMTM4ODgsNDIuNjkwODIwMyBsIDQ1LjA0Nzg3MzIsNDIuMDc2NzU3OCBsIDQ1LjA0Nzg3MzIsNDEuMjMyNDIxOSBsIDQ1LjIwMTM4ODgsNDAuNTQxNjAxNiBsIDQ1LjY2MTkzNTcsNDAuMTU3ODEyNSBsIDQ2LjA0NTcyNDcsMzkuNjk3MjY1NiBsIDQ2LjUwNjI3MTYsMzkuMzEzNDc2NiBsIDQ2Ljg5MDA2MDcsMzguODUyOTI5NyBsIDQ3LjUwNDEyMzIsMzguNjIyNjU2MyBsIDQ4LjE5NDk0MzUsMzguNDY5MTQwNiBsIDQ4LjU3ODczMjYsMzguMDA4NTkzOCBsIDQ5LjAzOTI3OTQsMzcuNjI0ODA0NyBsIDQ5LjQyMzA2ODUsMzcuMTY0MjU3OCBsIDQ5Ljg4MzYxNTQsMzYuNzgwNDY4OCBsIDUwLjI2NzQwNDQsMzYuMzE5OTIxOSBsIDUwLjcyNzk1MTMsMzUuOTM2MTMyOCBsIDUwLjk1ODIyNDcsMzUuMzIyMDcwMyBsIDUwLjk1ODIyNDcsMzQuNjk5ODY2MyBsIDUwLjk1ODIyNDcsMzMuOTE1MDQ0MSBsIDUwLjk1ODIyNDcsMzIuOTQwMDEwOSBsIDUwLjk1ODIyNDcsMzIuMzA3NTU0NCBsIDUwLjk1ODIyNDcsMzEuNDY4MzM1MiBsIDUwLjk1ODIyNDcsMzAuNTg5MTIzMiBsIDUwLjk1ODIyNDcsMjkuODEzODc0MiBsIDUwLjk1ODIyNDcsMjkuMDIzNjg0IGwgNTAuOTU4MjI0NywyOC4xMjM4MDQ3IGwgNTAuOTU4MjI0NywyNy4wMDM5MjExIGwgNTAuOTU4MjI0NywyNi4wMzQzNzUgYyA1MC45NTgyMjQ5LDI1Ljg0NjAxMjUgNTEuMTExNzQwNCwyNS4zNDM1NTQ3IDUxLjExMTc0MDQsMjUuMzQzNTU0NyBsIDUxLjU3MjI4NzIsMjQuOTU5NzY1NiBsIDUxLjgwMjU2MDcsMjQuMzQ1NzAzMSBsIDUxLjgwMjU2MDcsMjMuNTAxMzY3MiBsIDUxLjU3MjI4NzIsMjIuODEwNTQ2OSBsIDUxLjExMTc0MDQsMjIuNDI2NzU3OCBsIDUwLjk1ODIyNDcsMjEuODEyNjk1MyBsIDUwLjk1ODIyNDcsMjEuMDI5MDkxMiBsIDUwLjk1ODIyNDcsMjAuMzM1ODg1NCBsIDUwLjk1ODIyNDcsMTkuNjM0NzE2NyBsIDUwLjk1ODIyNDcsMTguOTU3NTI1OCBsIDUwLjk1ODIyNDcsMTguMDAwNTY1MyBsIDUwLjk1ODIyNDcsMTcuMzU1MzE0OCBsIDUwLjk1ODIyNDcsMTYuNzQ2Njc5NyBsIDUwLjcyNzk1MTMsMTYuMDU1ODU5NCBsIDUwLjI2NzQwNDQsMTUuNjcyMDcwMyBsIDUwLjExMzg4ODgsMTUuMDU4MDA3OCBsIDQ5Ljg4MzYxNTQsMTQuMzY3MTg3NSBsIDQ5LjQyMzA2ODUsMTMuOTgzMzk4NCBsIDQ5LjAzOTI3OTQsMTMuNTIyODUxNiBsIDQ4LjU3ODczMjYsMTMuMTM5MDYyNSBsIDQ4LjE5NDk0MzUsMTIuNjc4NTE1NiBsIDQ3LjczNDM5NjYsMTIuMjk0NzI2NiBsIDQ3LjM1MDYwNzYsMTEuODM0MTc5NyBsIDQ2Ljg5MDA2MDcsMTEuNDUwMzkwNiBsIDQ2LjUwNjI3MTYsMTAuOTg5ODQzOCBsIDQ1LjgxNTQ1MTMsMTAuNzU5NTcwMyBsIDQ1LjIwMTM4ODgsMTAuNjA2MDU0NyBsIDQ1LjIwMTM4ODgsMTAuMTQ1NTA3OCBsIDQ1LjY2MTkzNTcsOS43NjE3MTg3NSBsIDQ1Ljg5MjIwOTEsOS4xNDc2NTYyNSBsIDQ1Ljg5MjIwOTEsOC4zMDMzMjAzMSBsIDQ2LjA0NTcyNDcsNy42MTI1IGwgNDYuNTA2MjcxNiw3LjIyODcxMDk0IGwgNDYuODkwMDYwNyw2Ljc2ODE2NDA2IGwgNDcuMzUwNjA3Niw2LjM4NDM3NSBsIDQ3LjczNDM5NjYsNS45MjM4MjgxMyBsIDQ4LjM0ODQ1OTEsNS42OTM1NTQ2OSBsIDQ5LjAzOTI3OTQsNS41NDAwMzkwNiBsIDQ5LjQyMzA2ODUsNS4wNzk0OTIxOSBsIDUwLjAzNzEzMSw0Ljg0OTIxODc1IGwgNTAuNzI3OTUxMyw1LjA3OTQ5MjE5IGwgNTEuMTExNzQwNCw1LjU0MDAzOTA2IGwgNTEuNzI1ODAyOSw1LjY5MzU1NDY5IGwgNTIuNTcwMTM4OCw1LjY5MzU1NDY5IGwgNTMuMjYwOTU5MSw1LjkyMzgyODEzIGwgNTMuNjQ0NzQ4Miw2LjM4NDM3NSBsIDU0LjI1ODgxMDcsNi41Mzc4OTA2MyBsIDU0Ljk0OTYzMSw2Ljc2ODE2NDA2IGwgNTUuMzMzNDIwMSw3LjIyODcxMDk0IGwgNTUuNzkzOTY2OSw3LjYxMjUgbCA1Ni4xNzc3NTYsOC4wNzMwNDY4OCBsIDU2LjYzODMwMjksOC40NTY4MzU5NCBsIDU3LjAyMjA5MTksOC45MTczODI4MSBsIDU3LjQ4MjYzODgsOS4zMDExNzE4OCBsIDU3Ljg2NjQyNzksOS43NjE3MTg3NSBsIDU4LjMyNjk3NDcsMTAuMTQ1NTA3OCBsIDU4LjcxMDc2MzgsMTAuNjA2MDU0NyBsIDU5LjE3MTMxMDcsMTAuNjA2MDU0NyBsIDU5LjU1NTA5OTcsMTAuMTQ1NTA3OCBsIDYwLjE2OTE2MjIsOS45MTUyMzQzOCBsIDYxLjAxMzQ5ODIsOS45MTUyMzQzOCBsIDYxLjcwNDMxODUsOS43NjE3MTg3NSBsIDYyLjA4ODEwNzYsOS4zMDExNzE4OCBsIDYyLjU0ODY1NDQsOC45MTczODI4MSBsIDYyLjkzMjQ0MzUsOC40NTY4MzU5NCBsIDYzLjM5Mjk5MDQsOC4wNzMwNDY4OCBsIDYzLjc3Njc3OTQsNy42MTI1IGwgNjQuMzkwODQxOSw3LjM4MjIyNjU2IGwgNjUuMDgxNjYyMiw3LjIyODcxMDk0IGwgNjUuNDY1NDUxMyw2Ljc2ODE2NDA2IGwgNjYuMDc5NTEzOCw2LjUzNzg5MDYzIGwgNjYuNjQ1NjY3OSw2LjUzNzg5MDYzIGwgNjcuMzAyMTAyLDYuNTM3ODkwNjMgbCA2OC4wNTkxODg3LDYuNTM3ODkwNjMgbCA2OC42MTI1MjE2LDYuNTM3ODkwNjMgbCA2OS4zMDMzNDE5LDYuMzg0Mzc1IGwgNjkuNjg3MTMxLDUuOTIzODI4MTMgbCA3MC4zMDExOTM1LDUuNjkzNTU0NjkgbCA3MC45OTA4MTk3LDUuNjkzNTU0NjkgbCA3MS41NDgxMjIyLDUuNjkzNTU0NjkgbCA3Mi4wNDM2MDE1LDUuNjkzNTU0NjkgbCA3Mi43ODU4MzY0LDUuNjkzNTU0NjkgbCA3My42Nzg1MzcyLDUuNjkzNTU0NjkgbCA3NC4zNjkzNTc2LDUuOTIzODI4MTMgbCA3NC43NTMxNDY2LDYuMzg0Mzc1IGwgNzUuMzY3MjA5MSw2LjUzNzg5MDYzIGwgNzYuMDU4MDI5NCw2Ljc2ODE2NDA2IGwgNzYuNDQxODE4NSw3LjIyODcxMDk0IGwgNzcuMDU1ODgxLDcuMzgyMjI2NTYgbCA3Ny43NDY3MDEzLDcuNjEyNSBsIDc4LjEzMDQ5MDQsOC4wNzMwNDY4OCBsIDc4LjU5MTAzNzIsOC40NTY4MzU5NCBsIDc4Ljk3NDgyNjMsOC45MTczODI4MSBsIDc5LjQzNTM3MzIsOS4zMDExNzE4OCBsIDc5LjgxOTE2MjIsOS43NjE3MTg3NSBsIDgwLjI3OTcwOTEsMTAuMTQ1NTA3OCBsIDgwLjUwOTk4MjYsMTAuODM2MzI4MSBsIDgwLjY2MzQ5ODIsMTEuNDUwMzkwNiBsIDgxLjEyNDA0NTEsMTEuODM0MTc5NyBsIDgxLjM1NDMxODUsMTIuNTI1IGwgODEuNTA3ODM0MSwxMy4xMzkwNjI1IGwgODEuOTY4MzgxLDEzLjUyMjg1MTYgbCA4Mi4xOTg2NTQ0LDE0LjIxMzY3MTkgbCA4Mi4xOTg2NTQ0LDE0Ljk0MDA5ODkgbCA4Mi4xOTg2NTQ0LDE1Ljk3OTEwMjUgbCA4Mi4xOTg2NTQ0LDE2Ljc0NjY3OTcgbCA4MS45NjgzODEsMTcuMzYwNzQyMiBsIDgxLjUwNzgzNDEsMTcuNzQ0NTMxMyBsIDgxLjUwNzgzNDEsMTguMjA1MDc4MSBsIDgxLjk2ODM4MSwxOC41ODg4NjcyIGwgODEuOTY4MzgxLDE5LjA0OTQxNDEgbCA4MS41MDc4MzQxLDE5LjQzMzIwMzEgbCA4MS4zNTQzMTg1LDIwLjEyNDAyMzQgbCA4MS4zNTQzMTg1LDIwLjk2NTM4OTMgbCA4MS4zNTQzMTg1LDIxLjgxMjY5NTMgYyA4MS4zNTQzMTkzLDIyLjQyNjc1ODIgODEuNTA3ODM0MSwyMi40MjY3NTc4IDgxLjUwNzgzNDEsMjIuNDI2NzU3OCBsIDgxLjk2ODM4MSwyMi44MTA1NDY5IGwgODIuMTk4NjU0NCwyMy41MDEzNjcyIGwgODIuMTk4NjU0NCwyNC4zNDIzMDc2IGwgODIuMTk4NjU0NCwyNS4yMTkyODI4IGwgODIuMTk4NjU0NCwyNS45Mjk1NzcyIGwgODIuMTk4NjU0NCwyNy4wNzEzODA3IGwgODIuMTk4NjU0NCwyNy42ODQ2OTA4IGwgODIuMTk4NjU0NCwyOC41NjczODI4IGMgODIuMTk4NjUyOSwyOS4xODE0NDQ1IDgyLjM1MjE3MDEsMjkuMTgxNDQ1MyA4Mi4zNTIxNzAxLDI5LjE4MTQ0NTMgbCA4Mi44MTI3MTY5LDI5LjU2NTIzNDQgbCA4My4wNDI5OTA0LDMwLjI1NjA1NDcgbCA4My4wNDI5OTA0LDMxLjEwMDM5MDYgbCA4Mi44MTI3MTY5LDMxLjcxNDQ1MzEgbCA4Mi4zNTIxNzAxLDMyLjA5ODI0MjIgbCA4Mi4xOTg2NTQ0LDMyLjc4OTA2MjUgbCA4Mi4zNTIxNzAxLDMzLjQwMzEyNSBsIDgyLjgxMjcxNjksMzMuNzg2OTE0MSBsIDgzLjA0Mjk5MDQsMzQuNDc3NzM0NCBsIDgzLjE5NjUwNiwzNS4wOTE3OTY5IGwgODMuNjU3MDUyOSwzNS40NzU1ODU5IGwgODQuMDQwODQxOSwzNS45MzYxMzI4IGwgODQuNjU0OTA0NCwzNi4wODk2NDg0IGwgODUuMzQ1NzI0NywzNi4zMTk5MjE5IGwgODUuNzI5NTEzOCwzNi43ODA0Njg4IGwgODYuMzQzNTc2MywzNi45MzM5ODQ0IGwgODcuMDM0Mzk2NiwzNy4xNjQyNTc4IGwgODcuMjY0NjcwMSwzNy44NTUwNzgxIGwgODcuNDE4MTg1NywzOC40NjkxNDA2IGwgODcuODc4NzMyNiwzOC44NTI5Mjk3IGwgODguMTA5MDA2LDM5LjU0Mzc1IGwgODcuODc4NzMyNiw0MC4xNTc4MTI1IGwgODcuNDE4MTg1Nyw0MC41NDE2MDE2IHogbSAxMzAuNDg4OTI0LDMyLjA5ODI0MjIgbCAxMzAuMzM1NDA4LDMyLjc4OTA2MjUgbCAxMzAuMzM1NDA4LDMzLjY4ODY2MDYgbCAxMzAuMzM1NDA4LDM0LjQ3NzczNDQgbCAxMzAuMTA1MTM1LDM1LjA5MTc5NjkgbCAxMjkuNjQ0NTg4LDM1LjQ3NTU4NTkgbCAxMjkuNDkxMDcyLDM2LjE2NjQwNjMgbCAxMjkuNDkxMDcyLDM2Ljk3NDI2MzEgbCAxMjkuNDkxMDcyLDM3Ljg1NTA3ODEgbCAxMjkuMjYwNzk5LDM4LjQ2OTE0MDYgbCAxMjguODAwMjUyLDM4Ljg1MjkyOTcgbCAxMjguNjQ2NzM2LDM5LjU0Mzc1IGwgMTI4LjQxNjQ2Myw0MC4xNTc4MTI1IGwgMTI3Ljk1NTkxNiw0MC41NDE2MDE2IGwgMTI3LjU3MjEyNyw0MS4wMDIxNDg0IGwgMTI3LjExMTU4LDQxLjM4NTkzNzUgbCAxMjYuNzI3NzkxLDQxLjg0NjQ4NDQgbCAxMjYuMjY3MjQ0LDQyLjIzMDI3MzQgbCAxMjUuODgzNDU1LDQyLjY5MDgyMDMgbCAxMjUuNDIyOTA4LDQzLjA3NDYwOTQgbCAxMjUuMDM5MTE5LDQzLjUzNTE1NjMgbCAxMjQuNTc4NTcyLDQzLjkxODk0NTMgbCAxMjQuMTk0NzgzLDQ0LjM3OTQ5MjIgbCAxMjMuNTAzOTYzLDQ0LjUzMzAwNzggbCAxMjIuODg5OSw0NC43NjMyODEzIGwgMTIyLjUwNjExMSw0NS4yMjM4MjgxIGwgMTIxLjgxNTI5MSw0NS4zNzczNDM4IGwgMTIwLjk3MDk1NSw0NS4zNzczNDM4IGwgMTIwLjM1Njg5Miw0NS42MDc2MTcyIGwgMTE5Ljk3MzEwMyw0Ni4wNjgxNjQxIGwgMTE5LjI4MjI4Myw0Ni4yMjE2Nzk3IGwgMTE4LjYyNDE4NSw0Ni4yMjE2Nzk3IGwgMTE3Ljk0MTU2Niw0Ni4yMjE2Nzk3IGwgMTE3LjExODUwNyw0Ni4yMjE2Nzk3IGwgMTE2LjQ0MjA0Myw0Ni4yMjE2Nzk3IGwgMTE1LjcwNzU0OCw0Ni4yMjE2Nzk3IGwgMTE1LjA2MDYwMyw0Ni4yMjE2Nzk3IGwgMTE0LjQ0NjU0MSw0Ni4wNjgxNjQxIGwgMTE0LjA2Mjc1Miw0NS42MDc2MTcyIGwgMTEzLjM3MTkzMiw0NS4zNzczNDM4IGwgMTEyLjUyNzU5Niw0NS4zNzczNDM4IGwgMTExLjkxMzUzMyw0NS4yMjM4MjgxIGwgMTExLjUyOTc0NCw0NC43NjMyODEzIGwgMTEwLjgzODkyNCw0NC41MzMwMDc4IGwgMTEwLjIyNDg2MSw0NC4zNzk0OTIyIGwgMTA5Ljg0MTA3Miw0My45MTg5NDUzIGwgMTA5LjM4MDUyNSw0My41MzUxNTYzIGwgMTA4Ljk5NjczNiw0My4wNzQ2MDk0IGwgMTA4LjUzNjE4OSw0Mi42OTA4MjAzIGwgMTA4LjE1MjQsNDIuMjMwMjczNCBsIDEwNy42OTE4NTMsNDEuODQ2NDg0NCBsIDEwNy4zMDgwNjQsNDEuMzg1OTM3NSBsIDEwNi44NDc1MTcsNDEuMDAyMTQ4NCBsIDEwNi40NjM3MjgsNDAuNTQxNjAxNiBsIDEwNi4wMDMxODIsNDAuMTU3ODEyNSBsIDEwNS42MTkzOTIsMzkuNjk3MjY1NiBsIDEwNS4xNTg4NDYsMzkuMzEzNDc2NiBsIDEwNS4wMDUzMywzOC42OTk0MTQxIGwgMTA0Ljc3NTA1NywzOC4wMDg1OTM4IGwgMTA0LjMxNDUxLDM3LjYyNDgwNDcgbCAxMDQuMTYwOTk0LDM3LjAxMDc0MjIgbCAxMDQuMTYwOTk0LDM2LjE2NjQwNjMgbCAxMDMuOTMwNzIxLDM1LjQ3NTU4NTkgbCAxMDMuNDcwMTc0LDM1LjA5MTc5NjkgbCAxMDMuMzE2NjU4LDM0LjQ3NzczNDQgbCAxMDMuMzE2NjU4LDMzLjc0NDA1MzQgbCAxMDMuMzE2NjU4LDMyLjc4OTA2MjUgYyAxMDMuMzE2NjU4LDMyLjQzMTE3MDcgMTAzLjQ3MDE3NCwzMi4wOTgyNDIyIDEwMy40NzAxNzQsMzIuMDk4MjQyMiBsIDEwMy45MzA3MjEsMzEuNzE0NDUzMSBsIDEwNC4xNjA5OTQsMzEuMTAwMzkwNiBsIDEwNC4xNjA5OTQsMzAuMzk2MjYzOSBsIDEwNC4xNjA5OTQsMjkuNTE1MDk0MSBsIDEwNC4xNjA5OTQsMjguNDMwNjU3NCBsIDEwNC4xNjA5OTQsMjcuNzIzMDQ2OSBjIDEwNC4xNjA5OTQsMjcuNTQ0MjkxNyAxMDQuMzE0NTEsMjcuMDMyMjI2NiAxMDQuMzE0NTEsMjcuMDMyMjI2NiBsIDEwNC43NzUwNTcsMjYuNjQ4NDM3NSBsIDEwNS4wMDUzMywyNi4wMzQzNzUgbCAxMDQuNzc1MDU3LDI1LjM0MzU1NDcgbCAxMDQuMzE0NTEsMjQuOTU5NzY1NiBsIDEwNC4xNjA5OTQsMjQuMzQ1NzAzMSBsIDEwNC4zMTQ1MSwyMy42NTQ4ODI4IGwgMTA0Ljc3NTA1NywyMy4yNzEwOTM4IGwgMTA1LjAwNTMzLDIyLjY1NzAzMTMgbCAxMDUuMDA1MzMsMjEuNTQxMTYzNCBsIDEwNS4wMDUzMywyMC42MjcxNDE2IGwgMTA1LjAwNTMzLDE5Ljk5MjI0MjIgbCAxMDUuMDA1MzMsMTkuMjU2NzA4NSBsIDEwNS4wMDUzMywxOC40OTU1OTY1IGwgMTA1LjAwNTMzLDE3LjY0Nzc1ODQgbCAxMDUuMDA1MzMsMTYuODAwMDAwMiBsIDEwNS4wMDUzMywxNS45MDIzNDM4IGwgMTA0Ljc3NTA1NywxNS4yMTE1MjM0IGwgMTA0LjMxNDUxLDE0LjgyNzczNDQgbCAxMDQuMTYwOTk0LDE0LjIxMzY3MTkgbCAxMDMuOTMwNzIxLDEzLjUyMjg1MTYgbCAxMDMuNDcwMTc0LDEzLjEzOTA2MjUgbCAxMDMuMDg2Mzg1LDEyLjY3ODUxNTYgbCAxMDIuMzk1NTY0LDEyLjQ0ODI0MjIgbCAxMDEuNzgxNTAyLDEyLjI5NDcyNjYgbCAxMDEuMzk3NzEzLDExLjgzNDE3OTcgbCAxMDAuNzA2ODkyLDExLjYwMzkwNjMgbCAxMDAuMDUzMzU3LDExLjYwMzkwNjMgbCA5OS40Mjg4NDgzLDExLjYwMzkwNjMgbCA5OC44MjAwMDY3LDExLjYwMzkwNjMgbCA5OC4xNzM4ODQ2LDExLjYwMzkwNjMgbCA5Ny41NTk4MjIxLDExLjQ1MDM5MDYgbCA5Ny4xNzYwMzMxLDEwLjk4OTg0MzggbCA5Ni43MTU0ODYyLDEwLjYwNjA1NDcgbCA5Ni4zMzE2OTcxLDEwLjE0NTUwNzggbCA5NS44NzExNTAzLDkuNzYxNzE4NzUgbCA5NS40ODczNjEyLDkuMzAxMTcxODggbCA5NS4wMjY4MTQzLDguOTE3MzgyODEgbCA5NC44NzMyOTg3LDguMzAzMzIwMzEgbCA5NS4wMjY4MTQzLDcuNjEyNSBsIDk1LjQ4NzM2MTIsNy4yMjg3MTA5NCBsIDk1Ljg3MTE1MDMsNi43NjgxNjQwNiBsIDk2LjMzMTY5NzEsNi4zODQzNzUgbCA5Ni43MTU0ODYyLDUuOTIzODI4MTMgbCA5Ny4xNzYwMzMxLDUuNTQwMDM5MDYgbCA5Ny41NTk4MjIxLDUuMDc5NDkyMTkgbCA5OC4xNzM4ODQ2LDQuODQ5MjE4NzUgbCA5OC43MzM3NjAxLDQuODQ5MjE4NzUgbCA5OS4yMTgzODc1LDQuODQ5MjE4NzUgbCA5OS43MDE5NzU4LDQuODQ5MjE4NzUgbCAxMDAuMTc4NTMsNC44NDkyMTg3NSBsIDEwMC42MDQwODgsNC44NDkyMTg3NSBsIDEwMS40MDM4MDcsNC44NDkyMTg3NSBsIDEwMi4wMDI2NTcsNC44NDkyMTg3NSBsIDEwMi4zOTU1NjQsNC44NDkyMTg3NSBsIDEwMy4wODYzODUsNS4wNzk0OTIxOSBsIDEwMy40NzAxNzQsNS41NDAwMzkwNiBsIDEwMy45MzA3MjEsNS41NDAwMzkwNiBsIDEwNC4zMTQ1MSw1LjA3OTQ5MjE5IGwgMTA0Ljc3NTA1Nyw0LjY5NTcwMzEzIGwgMTA1LjAwNTMzLDQuMDgxNjQwNjMgbCAxMDUuMTU4ODQ2LDMuMzkwODIwMzEgbCAxMDUuNjE5MzkyLDMuMDA3MDMxMjUgbCAxMDUuODQ5NjY2LDIuMzkyOTY4NzUgbCAxMDYuMDAzMTgyLDEuNzAyMTQ4NDQgbCAxMDYuNjE3MjQ0LDEuNDcxODc1IGwgMTA3LjMwODA2NCwxLjMxODM1OTM4IGwgMTA3LjUzODMzOCwwLjcwNDI5Njg3NSBsIDEwNy42OTE4NTMsMC4wMTM0NzY1NjI1IGwgMTA4LjE1MjQsMC4wMTM0NzY1NjI1IGwgMTA4LjM4MjY3NCwwLjcwNDI5Njg3NSBsIDEwOC41MzYxODksMS4zMTgzNTkzOCBsIDEwOC45OTY3MzYsMS43MDIxNDg0NCBsIDEwOS4zODA1MjUsMi4xNjI2OTUzMSBsIDEwOS44NDEwNzIsMi41NDY0ODQzOCBsIDExMC4wNzEzNDYsMy4yMzczMDQ2OSBsIDExMC4wNzEzNDYsNC4wODE2NDA2MyBsIDExMC4yMjQ4NjEsNC42OTU3MDMxMyBsIDExMC42ODU0MDgsNS4wNzk0OTIxOSBsIDExMS4wNjkxOTcsNS41NDAwMzkwNiBsIDExMS42ODMyNiw1LjY5MzU1NDY5IGwgMTEyLjM3NDA4LDUuOTIzODI4MTMgbCAxMTIuNzU3ODY5LDYuMzg0Mzc1IGwgMTEzLjIxODQxNiw2LjM4NDM3NSBsIDExMy42MDIyMDUsNS45MjM4MjgxMyBsIDExNC4yMTYyNjcsNS42OTM1NTQ2OSBsIDExNS4wMDQ2MjcsNS42OTM1NTQ2OSBsIDExNS45MDQ5MzksNS42OTM1NTQ2OSBsIDExNi41OTU3Niw1LjkyMzgyODEzIGwgMTE2Ljk3OTU0OSw2LjM4NDM3NSBsIDExNy41OTM2MTEsNi41Mzc4OTA2MyBsIDExOC4yNTE3ODIsNi41Mzc4OTA2MyBsIDExOC43MzIwOTMsNi41Mzc4OTA2MyBsIDExOS4zODkwNTQsNi41Mzc4OTA2MyBjIDExOS42MTAwOTIsNi41Mzc4OTA2MyAxMTkuODMxMTMsNi41Mzc4OTA2MyAxMjAuMDUyMTY5LDYuNTM3ODkwNjMgYyAxMjAuMjAzNzY1LDYuNTM3ODkwNjMgMTIwLjM4MTgwMyw2LjUzNzg5MDYzIDEyMC41NzI0NDYsNi41Mzc4OTA2MyBjIDEyMC43NjQ5MjIsNi41Mzc4OTA2MyAxMjAuOTU3Mzk4LDYuNTM3ODkwNjMgMTIxLjE0OTg3NCw2LjUzNzg5MDYzIGMgMTIxLjM1MjEzNSw2LjUzNzg5MDYzIDEyMS42MzU2ODYsNi41Mzc4OTA2MyAxMjEuOTAyNTEzLDYuNTM3ODkwNjMgbCAxMjIuNjU5NjI3LDYuNTM3ODkwNjMgYyAxMjIuODk3NTcsNi41Mzc4OTA2MyAxMjMuMzUwNDQ3LDYuNzY4MTY0MDYgMTIzLjM1MDQ0Nyw2Ljc2ODE2NDA2IGwgMTIzLjczNDIzNiw3LjIyODcxMDk0IGwgMTI0LjE5NDc4Myw3LjYxMjUgbCAxMjQuNTc4NTcyLDguMDczMDQ2ODggbCAxMjUuMDM5MTE5LDguNDU2ODM1OTQgbCAxMjUuMjY5MzkyLDkuMTQ3NjU2MjUgbCAxMjUuMDM5MTE5LDkuNzYxNzE4NzUgbCAxMjQuNTc4NTcyLDEwLjE0NTUwNzggbCAxMjQuMTk0NzgzLDEwLjYwNjA1NDcgbCAxMjMuNzM0MjM2LDEwLjk4OTg0MzggbCAxMjMuMzUwNDQ3LDExLjQ1MDM5MDYgbCAxMjIuNjU5NjI3LDExLjYwMzkwNjMgbCAxMjEuOTg4MiwxMS42MDM5MDYzIGwgMTIxLjM0OTM4NCwxMS42MDM5MDYzIGwgMTIwLjc4MDk4OCwxMS42MDM5MDYzIGMgMTIwLjU5MzYwMSwxMS42MDM5MDYzIDEyMC40MDYyMTQsMTEuNjAzOTA2MyAxMjAuMjE4ODI3LDExLjYwMzkwNjMgbCAxMTkuNDQyNDQ4LDExLjYwMzkwNjMgbCAxMTguNzc0MTM3LDExLjYwMzkwNjMgbCAxMTguMTc2MjQ2LDExLjYwMzkwNjMgbCAxMTcuNTIwMjQ1LDExLjYwMzkwNjMgbCAxMTYuODM2NjY3LDExLjYwMzkwNjMgbCAxMTYuMjM4Nzc2LDExLjYwMzkwNjMgYyAxMTYuMDUxODYzLDExLjYwMzkwNjMgMTE1Ljg2ODY0OSwxMS42MDM5MDYzIDExNS42OTA5MjIsMTEuNjAzOTA2MyBjIDExNS40OTA2MjIsMTEuNjAzOTA2MyAxMTUuMjk3MjkyLDExLjYwMzkwNjMgMTE1LjExMzQ5NCwxMS42MDM5MDYzIGwgMTE0LjUzNjA2NiwxMS42MDM5MDYzIGMgMTE0LjM3NDkwNywxMS42MDM5MDYzIDExNC4yMjY0NDEsMTEuNjAzOTA2MyAxMTQuMDkzMzIyLDExLjYwMzkwNjMgYyAxMTMuNjQ1OTY3LDExLjYwMzkwNjMgMTEzLjM3MTkzMiwxMS42MDM5MDYzIDExMy4zNzE5MzIsMTEuNjAzOTA2MyBsIDExMi43NTc4NjksMTEuODM0MTc5NyBsIDExMi4zNzQwOCwxMi4yOTQ3MjY2IGwgMTExLjkxMzUzMywxMi42Nzg1MTU2IGwgMTExLjUyOTc0NCwxMy4xMzkwNjI1IGwgMTExLjA2OTE5NywxMy41MjI4NTE2IGwgMTEwLjY4NTQwOCwxMy45ODMzOTg0IGwgMTEwLjIyNDg2MSwxNC4zNjcxODc1IGwgMTEwLjA3MTM0NiwxNS4wNTgwMDc4IGwgMTEwLjA3MTM0NiwxNS45OTAyODk2IGwgMTEwLjA3MTM0NiwxNi43NDY2Nzk3IGMgMTEwLjA3MTM0NSwxNy4yMTg5MjI4IDExMC4yMjQ4NjEsMTcuMzYwNzQyMiAxMTAuMjI0ODYxLDE3LjM2MDc0MjIgbCAxMTAuNjg1NDA4LDE3Ljc0NDUzMTMgbCAxMTAuNjg1NDA4LDE4LjIwNTA3ODEgbCAxMTAuMjI0ODYxLDE4LjU4ODg2NzIgbCAxMTAuMDcxMzQ2LDE5LjI3OTY4NzUgbCAxMTAuMDcxMzQ2LDIwLjE3MzA0ODMgYyAxMTAuMDcxMzQ2LDIwLjM5NTA1NTIgMTEwLjA3MTM0NiwyMC42MzQ2NzY2IDExMC4wNzEzNDYsMjAuODg3NzE5OCBsIDExMC4wNzEzNDYsMjEuNjU2OTg0OCBsIDExMC4wNzEzNDYsMjIuMzkwNTIwMiBsIDExMC4wNzEzNDYsMjMuMTg3MzYxOCBsIDExMC4wNzEzNDYsMjMuOTcwOTM0NiBsIDExMC4wNzEzNDYsMjQuNjUzNTUzNCBsIDExMC4wNzEzNDYsMjUuNDczNzM1MSBsIDExMC4wNzEzNDYsMjYuMjY1MzgxIGwgMTEwLjA3MTM0NiwyNi45MTMyMjk0IGwgMTEwLjA3MTM0NiwyNy41NDQyOTIxIGMgMTEwLjA3MTM0NiwyOC4xNzQyNTQyIDExMC4wNzEzNDYsMjguNTY3MzgyOCAxMTAuMDcxMzQ2LDI4LjU2NzM4MjggbCAxMTAuMjI0ODYxLDI5LjE4MTQ0NTMgbCAxMTAuNjg1NDA4LDI5LjU2NTIzNDQgbCAxMTAuOTE1NjgyLDMwLjI1NjA1NDcgbCAxMTAuOTE1NjgyLDMxLjEwMDM5MDYgbCAxMTEuMDY5MTk3LDMxLjcxNDQ1MzEgbCAxMTEuNTI5NzQ0LDMyLjA5ODI0MjIgbCAxMTEuNzYwMDE3LDMyLjc4OTA2MjUgbCAxMTEuNzYwMDE3LDMzLjYzMzM5ODQgbCAxMTEuOTEzNTMzLDM0LjI0NzQ2MDkgbCAxMTIuMzc0MDgsMzQuNjMxMjUgbCAxMTIuNjA0MzUzLDM1LjMyMjA3MDMgbCAxMTIuNzU3ODY5LDM1LjkzNjEzMjggbCAxMTMuMzcxOTMyLDM2LjA4OTY0ODQgbCAxMTQuMDYyNzUyLDM2LjMxOTkyMTkgbCAxMTQuNDQ2NTQxLDM2Ljc4MDQ2ODggbCAxMTUuMDYwNjAzLDM2LjkzMzk4NDQgbCAxMTUuNzQyMzE4LDM2LjkzMzk4NDQgbCAxMTYuMjI2ODY2LDM2LjkzMzk4NDQgbCAxMTYuNjgzOTE3LDM2LjkzMzk4NDQgbCAxMTcuMDQxODUyLDM2LjkzMzk4NDQgbCAxMTcuNzI3MzQ4LDM2LjkzMzk4NDQgbCAxMTguNDM3OTQ3LDM2LjkzMzk4NDQgbCAxMTkuMTI4NzY3LDM2Ljc4MDQ2ODggbCAxMTkuNTEyNTU3LDM2LjMxOTkyMTkgbCAxMTkuOTczMTAzLDM1LjkzNjEzMjggbCAxMjAuMzU2ODkyLDM1LjQ3NTU4NTkgbCAxMjAuODE3NDM5LDM1LjA5MTc5NjkgbCAxMjEuMDQ3NzEzLDM0LjQ3NzczNDQgbCAxMjEuMjAxMjI4LDMzLjc4NjkxNDEgbCAxMjEuNjYxNzc1LDMzLjQwMzEyNSBsIDEyMS44OTIwNDksMzIuNzg5MDYyNSBsIDEyMS44OTIwNDksMzEuOTQ0NzI2NiBsIDEyMi4wNDU1NjQsMzEuMjUzOTA2MyBsIDEyMi41MDYxMTEsMzAuODcwMTE3MiBsIDEyMi43MzYzODUsMzAuMjU2MDU0NyBsIDEyMi44ODk5LDI5LjU2NTIzNDQgbCAxMjMuMzUwNDQ3LDI5LjE4MTQ0NTMgbCAxMjMuNzM0MjM2LDI4LjcyMDg5ODQgbCAxMjQuMzQ4Mjk5LDI4LjQ5MDYyNSBsIDEyNS4wMzkxMTksMjguMzM3MTA5NCBsIDEyNS40MjI5MDgsMjcuODc2NTYyNSBsIDEyNi4wMzY5NzEsMjcuNjQ2Mjg5MSBsIDEyNi44ODEzMDcsMjcuNjQ2Mjg5MSBsIDEyNy41NzIxMjcsMjcuODc2NTYyNSBsIDEyNy45NTU5MTYsMjguMzM3MTA5NCBsIDEyOC41Njk5NzgsMjguNDkwNjI1IGwgMTI5LjI2MDc5OSwyOC43MjA4OTg0IGwgMTI5LjY0NDU4OCwyOS4xODE0NDUzIGwgMTMwLjEwNTEzNSwyOS41NjUyMzQ0IGwgMTMwLjMzNTQwOCwzMC4yNTYwNTQ3IGwgMTMwLjQ4ODkyNCwzMC44NzAxMTcyIGwgMTMwLjk0OTQ3MSwzMS4yNTM5MDYzIGwgMTMwLjk0OTQ3MSwzMS43MTQ0NTMxIGwgMTMwLjQ4ODkyNCwzMi4wOTgyNDIyIHogXCI7XG52YXIgbGV0dHJlUyA9IFwibSAxMzguODExNTE3LDM1LjEwNjEwMTEgYyAxMzguNjM3ODEsMzUuNzM2NTExMSAxMzguNjA2NTY0LDM1LjgwMTA2MTQgMTM4LjE5ODQ1NSwzNi4yMDkxNzA0IGMgMTM4LjA1Njc1MywzNi4zNTA4NzIyIDEzNy45OTkxNTUsMzcuMDIxNTkyOCAxMzcuOTQ4Nzk4LDM3LjIyMzAxOTkgYyAxMzcuODY5NTM3LDM3LjU0MDA2MjIgMTM3LjkyNjQxMSwzNy44MjkyMzMxIDEzNy45NDg3OTgsMzguMDc3NTQ3NSBjIDEzNy45ODMwNzgsMzguNDU3NzgwNiAxMzcuOTI3NzY2LDM4Ljk0MTUzNjcgMTM4LjEyNjAzNywzOS4yMjE3NDc5IGMgMTM4LjI0MDcxMywzOS4zODM4MTY0IDEzOC40Mzg1OTYsMzkuNDMzODM1IDEzOC41NjA1NDQsMzkuNjIwMTk2NCBjIDEzOC43NDcwNjYsMzkuOTA1MjQwOCAxMzguNzk1OTEsNDAuNDU4NDUyNyAxMzguOTUxNTk5LDQwLjgzMzkzMyBjIDEzOS4wNjExMzEsNDEuMDk4MDkyMiAxMzkuMjYzOTE4LDQxLjIwNDU5NTggMTM5LjQyOTU1Niw0MS4zNzAyMzQzIGMgMTM5LjU5MjY3Miw0MS41MzMzNDk1IDEzOS42MjY0MTIsNDEuOTM0NTM2MyAxMzkuNjIyMTY5LDQyLjQyNDU4NzYgYyAxMzkuNjE4MzIsNDIuODY5MDYyNSAxMzkuNTgzMjI2LDQzLjM4NjY0MSAxMzkuNTg1MDc2LDQzLjg2NTk5NDcgYyAxMzkuNTg2NTc5LDQ0LjI1NTQ1OTQgMTM5LjYxMjQ3MSw0NC42MTk2OTEyIDEzOS42OTkzMjQsNDQuODk4OTgwNiBjIDEzOS43OTg1MTgsNDUuMjE3OTU1NyAxNDAuMDU4NTMzLDQ1LjMzNjg1MTUgMTQwLjI0OTM4Miw0NS41MTYyNDM4IGMgMTQwLjQ2ODM3Miw0NS43MjIwODc3IDE0MC42Njc2OTgsNDUuOTUyNzQwNSAxNDAuODgzNzE0LDQ2LjIwNjQxOTYgYyAxNDEuMjg3OTg3LDQ2LjY4MTE3OTEgMTQxLjM2MzU2Myw0Ni44MDQyNjYyIDE0MS43NTg3MzMsNDYuODg4OTQ1MiBjIDE0Mi4xNTM5MDMsNDYuOTczNjI0MSAxNDIuMzg5NTcxLDQ2Ljg4ODk0NTIgMTQyLjU3NzMsNDYuODg4OTQ1MiBjIDE0Mi45MTkwMTgsNDYuODg4OTQ1MiAxNDMuMDYyNTg2LDQ2LjU0NDc5MjQgMTQzLjE3MDA1NSw0Ni40MzczMjI5IGMgMTQzLjUyNDQsNDYuMDgyOTc3NiAxNDQuMzAyNjEsNDYuMDg0OTMyNyAxNDQuNDU2MTA3LDQ1Ljk1NzY5MjggYyAxNDQuNjA5NjAzLDQ1LjgzMDQ1MjggMTQ0LjUzODI1Miw0NS43ODg4OTggMTQ0LjY4MDE2OSw0NS42NDY5ODExIGMgMTQ0Ljg2NzM5MSw0NS40NTk3NTkyIDE0NS4zNDA4LDQ1LjM0OTIzOTUgMTQ1LjkyMTA0Nyw0NS4yOTE5MzI5IGMgMTQ2LjQxNjg3OSw0NS4yNDI5NjMzIDE0Ni45OTA3MjMsNDUuMjMyODUgMTQ3LjUzMDY3Myw0NS4yNDY5MzY1IGMgMTQ4LjIzNDU5OSw0NS4yNjUzMDA4IDE0OC44ODA5MTcsNDUuMzI0Nzk1MiAxNDkuMjIxNjYxLDQ1LjM5Mjk0NDIgYyAxNDkuNzkzNTI3LDQ1LjUwNzMxNzIgMTQ5LjU5MjI5Myw0NS45NDA2Nzk0IDE1MC4zMTEzNDQsNDYuMDg0NDg5NiBjIDE1MC41OTA2MTUsNDYuMTQwMzQzNyAxNTAuOTgzNTA4LDQ2LjEzMjg1OTYgMTUxLjM4NTcwNSw0Ni4xMzM5ODMzIGMgMTUxLjc0MDQ2MSw0Ni4xMzQ5NzQ1IDE1Mi4xMDI0NTUsNDYuMTQyNjYyNiAxNTIuNDAwMTAxLDQ2LjIwNjQxOTYgYyAxNTMuMDM1MTk1LDQ2LjM0MjQ1OTYgMTUyLjg0MzI1OSw0Ni42OTg5NjQ5IDE1My4zMDM3OTcsNDYuODU1MzU2MyBjIDE1My43NDAxMDksNDcuMDAzNTIwOCAxNTQuMTI2NjU5LDQ3LjA2MDcxNzcgMTU0LjQ2Mjc2Nyw0Ny4wNTI1MzA3IGMgMTU0Ljg0ODg1OCw0Ny4wNDMxMjYzIDE1NS4xNjgzODksNDYuOTQ3NDQ1MSAxNTUuNDIwMzI5LDQ2LjgwNDI2NjIgYyAxNTUuNzAxMTUzLDQ2LjY0NDY3MjMgMTU1Ljg1NjU4Miw0Ni4zMDE4ODUyIDE1Ni4yMzg0NDQsNDYuMjA2NDE5NiBjIDE1Ni42Mzc1NSw0Ni4xMDY2NDMgMTU3LjA2MjQ2Myw0Ni4wNzczODY0IDE1Ny41MDExOTMsNDYuMDc5OTAwMyBjIDE1Ny45NjI1MTQsNDYuMDgyNTQzNiAxNTguNDM5MTExLDQ2LjEyMDMxMzMgMTU4LjkxNzA0Nyw0Ni4xNDgxNjA5IGMgMTU5LjM5NTE1OCw0Ni4xNzYwMTg3IDE1OS44NzQ2MDksNDYuMTkzOTQ3IDE2MC4zNDE0NDUsNDYuMTU2ODQ3NyBjIDE2MC42OTU2MjQsNDYuMTI4NzAxMSAxNjEuMDQyNTQzLDQ2LjA2ODg4MDkgMTYxLjM3NjEwNyw0NS45NTc2OTI4IGMgMTYxLjU1NjE2NSw0NS44OTc2NzM0IDE2MS41OTU1MDEsNDUuNjU0MjE1OSAxNjEuODU1OTU3LDQ1LjUxNjI0MzggYyAxNjIuMjE0MzU3LDQ1LjMyNjM4NzIgMTYyLjc2MDQ1LDQ1LjIwNDU3MjIgMTYyLjk5OTEyNyw0NS4xMjUwMTMzIGMgMTYzLjMwMTQ2OSw0NS4wMjQyMzI2IDE2My4yNTQyMTEsNDQuODE3Njg0MyAxNjMuNTQyNzk4LDQ0LjY3MzM5MSBjIDE2My44MTE5MjEsNDQuNTM4ODI5MyAxNjQuMTk4NzUxLDQ0LjUwNDAyOTQgMTY0LjUzNzQ2OSw0NC4zNjI4OTkxIGMgMTY0Ljg3NjE4Nyw0NC4yMjE3Njg3IDE2NS4yNDU2MjgsNDMuNzk1ODcxOCAxNjUuNTk1OTYsNDMuNDQ1NTQwNyBjIDE2NS42OTQwNzQsNDMuMzQ3NDI2MyAxNjYuMDg5NDI3LDQyLjk1MTU3NyAxNjYuNDE0MDMxLDQyLjU4NDYzMzcgYyAxNjYuNzM4NjM1LDQyLjIxNzY5MDMgMTY2LjY0OTM5OCw0MS42OTU4MzgxIDE2Ni44MjM4MSw0MS4zNzAyMzQzIGMgMTY2Ljk0MjU4Nyw0MS4xNDg0OTI4IDE2Ny4yMzE3MjgsNDEuMDcxNzIwMiAxNjcuMzQ2MTY0LDQwLjgzMzkzMyBjIDE2Ny41MTQ3OTcsNDAuNDgzNTI4NiAxNjcuNDcwNTY0LDM5Ljk1ODQzMzkgMTY3LjY0MjU0MiwzOS42MjAxOTY0IGMgMTY3Ljc1MjMzMiwzOS40MDQyNjgxIDE2OC4wNjQ3NTYsMzkuMzMwNDY0NiAxNjguMTY0NzMyLDM5LjEyNjIzMjggYyAxNjguMzM3NTc2LDM4Ljc3MzE0MzIgMTY4LjM4MjIzNSwzOC4yMTkyNDI1IDE2OC41MTczOTUsMzcuOTU0ODM3NSBjIDE2OC42NjM4MDIsMzcuNjY4NDI5NCAxNjguOTQ2MzUzLDM3LjUzNjY5NTEgMTY5LjAxMTM1OCwzNy40MDY2ODU3IGMgMTY5LjE2MTk5LDM3LjEwNTQyMTQgMTY5LjE5NjQxOCwzNi41ODM2MTk0IDE2OS4yMTYwMDEsMzYuMDUzNDMwMSBjIDE2OS4yMzY4MywzNS40ODk1MDU1IDE2OS4yNDA4NjYsMzQuOTE2MDkyNCAxNjkuMzUwMDc0LDM0LjU4ODQ2OTEgYyAxNjkuNDY5Nzc0LDM0LjIyOTM2OTMgMTY5LjgxNTgxLDM0LjI5MDY5MDIgMTY5Ljk3MTA1OCwzMy43Njg1MDI3IGMgMTcwLjEyNjMwNiwzMy4yNDYzMTUyIDE2OS42Nzg2MTMsMzIuNzY3NjM1NCAxNjkuNDA5MjksMzIuNDk4MzEyOSBjIDE2OS4yMDkzMDUsMzIuMjk4MzI3MyAxNjkuMzQxMDM1LDMxLjQ0Mzc2NTMgMTY4LjU3Mzg1LDMwLjg0NzA2NTkgYyAxNjguMzE5ODExLDMwLjY0OTQ4MDUgMTY4LjQ1NjU1NSwzMC4wMDg2NjY3IDE2OC4xNjQ3MzIsMjkuNjA1MzIzNCBjIDE2Ny45OTM0MzcsMjkuMzY4NTY5MSAxNjcuNzI0MTY5LDI5LjI4NTY3MzkgMTY3LjY0MjU0MiwyOS4wNDA3OTQ2IGMgMTY3LjQ1NTg5NiwyOC40ODA4NTY4IDE2Ny42MDkyNzUsMjguNDI3OTYyMiAxNjcuMzQ2MTY0LDI3Ljk4MjMwNCBjIDE2Ny4xMjQxOTQsMjcuNjA2MzI5MyAxNjYuNjIyNzgsMjcuMTQ1OTc1OSAxNjYuMjg3NjczLDI2Ljg1NTQwODggYyAxNjUuOTQ0NjM1LDI2LjU1Nzk2NDMgMTY1LjcxMTA1NCwyNi4yNDg0OTQ3IDE2NS40NTQzMzEsMjYuMDA2NDUzMSBjIDE2NS4xNDQ2MjgsMjUuNzE0NDYxMyAxNjUuMDQ1MDYxLDI1LjQ4NDU5NCAxNjQuNjc4NzY2LDI1LjI3MTEwMzQgYyAxNjQuMzQyMjcsMjUuMDc0OTgwNSAxNjMuOTIyMzk4LDI1LjA0NjY2MTcgMTYzLjU0Mjc5OCwyNC44NzczOTU1IGMgMTYzLjE5NDQ3OCwyNC43MjIwNzczIDE2My4wODM5NzIsMjQuMzU1MjA2MiAxNjIuNTA3NTM2LDI0LjIyODE4NzcgYyAxNjEuOTMxMSwyNC4xMDExNjkzIDE2MS40ODY1MDUsMjQuMjI3ODkzNiAxNjEuMDIzNDQyLDI0LjA3Mjk0MTggYyAxNjAuNjgzOTI2LDIzLjk1OTMzMTMgMTYwLjYxMzcyMiwyMy41NjAwMTY5IDE2MC4yNDcyMTQsMjMuNDM3ODQ3OCBjIDE1OS43ODgwMzIsMjMuMjg0Nzg3IDE1OS4xNDI4NzksMjMuMzA3OTAzMyAxNTguNDcwMTIyLDIzLjM0NzkxMTUgYyAxNTcuOTQ0NjEyLDIzLjM3OTE2MyAxNTcuNDAyMjYsMjMuNDIwNzIxMyAxNTYuOTE4NTQ1LDIzLjM5NjY2ODMgYyAxNTYuNTQwNTA3LDIzLjM3Nzg3MDEgMTU2LjE5ODI4NCwyMy4zMTg5OTcgMTU1LjkyNzkwOCwyMy4xODM4MDkxIGMgMTU1LjcwNjY1NiwyMy4wNzMxODI5IDE1NS40OTEwNjEsMjIuNjI1NDMwOCAxNTQuOTEyNDE4LDIyLjU0ODcxMzMgYyAxNTQuMzMzNzc2LDIyLjQ3MTk5NTcgMTUzLjgxNTI2OSwyMi41NDA4ODk5IDE1My40Nzc3NDMsMjIuMzk1MzczNCBjIDE1My4xNjI4NTcsMjIuMjU5NjE3NyAxNTMuMDU5NzYzLDIxLjg2Njk2MTIgMTUyLjkwNzY4LDIxLjgyODk0MDMgYyAxNTIuNjE5ODg1LDIxLjc1Njk5MTUgMTUyLjM0NDI5NiwyMS43MDg5ODM3IDE1Mi4wNzk3MTEsMjEuNjc3NzY3MSBjIDE1MS43MjM0ODksMjEuNjM1NzM4NiAxNTEuMzg3MjEzLDIxLjYyNDE0NjcgMTUxLjA2Nzk1MSwyMS42MjU1NDIyIGMgMTUwLjUzOTcxMiwyMS42Mjc4NTEyIDE1MC4wNTgwNSwyMS42NjU3MTQzIDE0OS42MDk2NzUsMjEuNjYwMDk2MSBjIDE0OS4xNjQyNTYsMjEuNjU0NTE0OSAxNDguNzUxNjg4LDIxLjYwNjAyNCAxNDguMzU4OTQyLDIxLjQzNzE0MTMgYyAxNDcuOTc1OTg3LDIxLjI3MjQ2OSAxNDcuODk4ODI2LDIxLjA5NTA1MjcgMTQ3LjIwODY1LDIwLjQxNzYxODMgYyAxNDYuOTk2NDc2LDIwLjIwOTM2MTQgMTQ2Ljc5OTE3NSwyMC4wMjY0NDY2IDE0Ni42MjA0ODIsMTkuODYyMzEzMyBjIDE0Ni4yMTc5MDksMTkuNDkyNTQxMyAxNDUuOTA5Nzg0LDE5LjIxODA5NDQgMTQ1LjczODgyNiwxOC45NjM5NTY4IGMgMTQ1LjQ5MTk4NCwxOC41OTcwMTM0IDE0NS41ODM1ODEsMTguMTQ1MzkwMiAxNDUuMzY4MTgyLDE3Ljc5MjU1OTcgYyAxNDUuMTUyNzgzLDE3LjQzOTcyOTIgMTQ1LjE2NTkyNiwxNy43OTg0OTUzIDE0NC44OTYwOTQsMTYuOTg4OTk5NSBjIDE0NC44MTA1OTMsMTYuNzMyNDk1MiAxNDUuMjg3NTk5LDE2LjYwMTQ0MzQgMTQ1LjM2ODE4MiwxNi4yOTg4MjM3IGMgMTQ1LjYzOTM3MiwxNS4yODA0MDY3IDE0NS42NDM5NTcsMTUuMTg1MDk2NiAxNDYuMTY1OTgyLDE0LjY1OTQyNDYgYyAxNDYuMzYwMTI3LDE0LjQ2MzkyMjggMTQ2LjQ5MDYyOSwxMy40OTc2NzY1IDE0Ni43MjYwOSwxMy4yNjIyMTY0IGMgMTQ2Ljk1ODczOSwxMy4wMjk1Njc0IDE0Ny4yMjI0MzksMTIuODIzMDkyMSAxNDcuNDU5MjY5LDEyLjU5MDU1NDMgYyAxNDcuNjkxOTA4LDEyLjM2MjEzMTQgMTQ3Ljg5ODYyMSwxMi4xMDg1NjAxIDE0OC4wMjQ1MDYsMTEuNzgwMzI4MyBjIDE0OC4yNzg1NDUsMTEuMTE3OTUxNSAxNDcuOTgyMTY3LDExLjA3OTEwNTMgMTQ4LjM1ODk0MiwxMC45MDUzMDg5IGMgMTQ4LjczNTcxNywxMC43MzE1MTI1IDE0OC43NzEzMDksMTEuMTYwOTQyIDE0OS4yMjE2NjEsMTEuMzk5MjcxNiBjIDE0OS41MDQ1NjcsMTEuNTQ4OTg3MyAxNDkuODIyMzE0LDExLjUyMDM2ODcgMTUwLjE5Nzk0MSwxMS41NTQ1MTY2IGMgMTUwLjU0NjQxLDExLjU4NjE5NTYgMTUxLjEyOTI5NywxMS41Nzc4Nzg3IDE1MS43MTYzOTcsMTEuNTY3NjI4MSBjIDE1Mi42MDYwMjQsMTEuNTUyMDk1NiAxNTMuNTA1MzIyLDExLjUzMjEyMzQgMTUzLjYxMzM0LDExLjY0MDE0MDkgYyAxNTMuNzY2NTA4LDExLjc5MzMwODUgMTUzLjk3MTU3NywxMi4wODY0NjU0IDE1NC4yMDYwOTYsMTIuMjAzNzI0OSBjIDE1NC42OTgzMTUsMTIuNDQ5ODM0NCAxNTUuMzA2NzMsMTIuNDI2NjYwOCAxNTUuOTIwMDczLDEyLjM4NTU5MzMgYyAxNTYuNjE0NzU1LDEyLjMzOTA3OTcgMTU3LjMxNTc1OSwxMi4yNjk2MTE1IDE1Ny44NjE0MiwxMi41NDI0NDIxIGMgMTU3Ljk4NDcyOCwxMi42MDQwOTYyIDE1OC4yMzMzNzMsMTIuOTU4MzQ5MyAxNTguMzQxNDY4LDEzLjAxMjM5NjkgYyAxNTkuMDI0Mjg0LDEzLjM1MzgwNDggMTU5LjI3ODIwNCwxMy4yNjE1Nzg4IDE1OS42Mzk2ODUsMTMuNDU4ODc1OCBjIDE1OS45NDIzMDgsMTMuNjI0MDQ3MyAxNjAuMTEzNTU5LDEzLjg5NjUzNjQgMTYwLjQ4NjQ3OCwxNC4yNzEwOTI3IGMgMTYwLjY4NjgwMSwxNC40NzIyOTQ2IDE2MC44ODY2NTIsMTQuNjg1OTgxIDE2MS4xODQ3LDE0Ljk4NDAyODUgYyAxNjEuMzY1NzgzLDE1LjE2NTExMTMgMTYxLjY5NTc2NSwxNS41MzcwNzY4IDE2Mi4wNTg5NTcsMTUuOTMyNjI4NyBjIDE2Mi4yNjczOSwxNi4xNTk2MzM2IDE2Mi40ODY3NjEsMTYuMzk0NDA2NyAxNjIuNjk1MjAyLDE2LjYwNTMyNjcgYyAxNjIuOTgxNTMsMTYuODk1MDYwOCAxNjMuMjQ3MjM0LDE3LjEzOTc4NDggMTYzLjQzNTYzMywxNy4yNTc1MzM4IGMgMTY0LjAwMDE2MiwxNy42MTAzNjQzIDE2NC4wMjEyODIsMTcuMzI1MjA1IDE2NC4zMzkzODcsMTcuNDExNTAyNiBjIDE2NC44MTAyMzgsMTcuNTM5MjM3OCAxNjUuMjY0MDQ3LDE3LjQ0MTQ2MjMgMTY1LjY1NzEyMywxNy4xOTg0ODA1IGMgMTY1Ljk5NTY3LDE2Ljk4OTIwNjMgMTY2LjI4OTE2NiwxNi42NzIyMTg5IDE2Ni41MDk2OTcsMTYuMjk4ODIzNyBjIDE2Ni43MDE5NTgsMTUuOTczMjk1IDE2Ni43NDEwMywxNS40ODgwNTUxIDE2Ni43MzQxMDQsMTQuOTEzOTkyNyBjIDE2Ni43MjgwODQsMTQuNDE0OTkzOCAxNjYuNjg3MzEsMTMuODQ4ODgyMiAxNjYuNjgyMTgyLDEzLjI2MjIxNjQgYyAxNjYuNjc5MDM2LDEyLjkwMjM4NDggMTY2LjY5NTE5NywxMi41MDA2ODk4IDE2Ni42ODc2OTUsMTIuMTI5NDMyNSBjIDE2Ni42Nzc0MzIsMTEuNjIxNTQwNSAxNjYuNjIyODgyLDExLjE3MDYxMyAxNjYuNDE0MDMxLDEwLjk2MTc2MTcgYyAxNjYuMTI4MDI4LDEwLjY3NTc1OTIgMTY2LjE1OTA3NiwxMC43OTA1NzA5IDE2NS45NjI0MDcsMTAuMzk3MjMyOSBjIDE2NS44NTYyMTIsMTAuMTg0ODQ0NCAxNjUuOTI0ODc3LDkuNTkyNDIyMDYgMTY1LjQ1NDMzMSw5LjE2OTM4Mjk3IGMgMTY0Ljk4Mzc4NSw4Ljc0NjM0Mzg4IDE2NC45NzQ0ODEsOC4zNjE4MDc5NyAxNjQuNDUyMjkyLDguMjY2MTM2ODUgYyAxNjQuMDY2MzA3LDguMTk1NDIwMDIgMTYzLjY4NTIyMyw4LjE0OTA5ODUyIDE2My4zMDYzOTcsOC4xMjAxMTIwOSBjIDE2Mi43OTE5NzMsOC4wODA3NTAyMyAxNjIuMjgxNzExLDguMDczMzU0MzUgMTYxLjc2ODk5NCw4LjA4MDI0NTI4IGMgMTYxLjI5NjAyNiw4LjA4NjYwMTk2IDE2MC44MjA5NjksOC4xMDUxMTYwNyAxNjAuMzM4NjI2LDguMTIxOTA5NzcgYyAxNTkuODUxNDE3LDguMTM4ODcyOTEgMTU5LjM1Njc3NSw4LjE1NDA4MDc2IDE1OC44NDkzNDUsOC4xNTMyMzExNiBjIDE1OC4zMDQ3MSw4LjE1MjMxOTI3IDE1Ny41MjM1LDguMTk3OTI0NzUgMTU3LjAyODc0LDguMDc0MjM0ODggYyAxNTYuNDc4NjczLDcuOTM2NzE3OTQgMTU2LjQxMDc1Nyw3LjUzNzk1MzUyIDE1NS45Mjc5MDgsNy4zNzcwMDQxIGMgMTU1LjQ0ODYyNCw3LjIxNzI0MjU1IDE1NC42NTk5MDksNy4zNzcwMDQxIDE1My43MDc4LDcuMjc4ODczMTMgYyAxNTMuMjQxMTY5LDcuMjMwNzc4ODUgMTUzLjE3NzQyNCw2LjkxMzY1MTk0IDE1Mi45MDc2OCw2LjcyODQ1NzYgYyAxNTIuNjM3Njg4LDYuNTQzMDkzMTQgMTUyLjIzMDU1LDYuNDE3OTQwODYgMTUxLjU2NjkyNCw2LjQ3NDQxOTY2IGMgMTUxLjEyMDc4LDYuNTEyMzg5MzMgMTUwLjgzMjA3LDYuNTczMjEyMTggMTUwLjM3MTk1Myw3LjEyMzYyNzcxIGMgMTQ5LjgwMjc3MSw3LjM2MzU1MjQxIDE0OS4yNzk2OTIsNy4zNzEwNTQzIDE0OS4wMTI0MzEsNy41MDQ2ODQ2MSBjIDE0OC44NzkwODksNy45NDkxNjE1OCAxNDguNDM0NzQ2LDguMDkzMzQ4MTYgMTQ3Ljk1MzAxOCw4LjE1MDcwMTkzIGMgMTQ3LjM2NDAyMyw4LjIyMDgyNjYzIDE0Ni43MTkxMzksOC4xNjExNDMyMiAxNDYuNTE4NDc1LDguMzYxODA3OTYgYyAxNDYuMzk0ODIxLDguNDg1NDYxNjkgMTQ2LjMxNDk1Niw4LjY1MzgyMjI2IDE0Ni4xNjU5ODIsOC44MDI3OTY1OSBjIDE0NS45NDg2NzgsOS4wMjAxMDA4NyAxNDUuNDQwMTkxLDguODkwNzkxNTggMTQ1LjEzMDMwOCw5LjA0NTczMjc2IGMgMTQ0Ljg4MDQzLDkuMTcwNjcxOTQgMTQ0LjY2MDI2LDkuNDQ4NjM3NDYgMTQ0LjQ1NjEwNyw5LjY1Mjc5MDcyIGMgMTQ0LjM0Njg0Myw5Ljc2MjA1NDEyIDE0My43OTgxNzIsOS42NzA4NDQ5MiAxNDMuMzA0MzM1LDEwLjAyOTczMTYgYyAxNDIuODEwNDk4LDEwLjM4ODYxODQgMTQyLjE1Nzg2NSwxMS4wNDMyNTEyIDE0Mi4wODk4NSwxMS4xMzU1NTgyIGMgMTQxLjk1NzYyOCwxMS4zMTUwMDExIDE0MS42NDM0MjQsMTEuNTEwODc0MyAxNDEuNTQzOSwxMS43MDk5MjI0IGMgMTQxLjEwMjc1MiwxMi41OTIyMTc2IDE0MS4zNTc4MjIsMTIuNDcwMzgwOSAxNDEuMjM2NDkzLDEyLjY5ODAxNDQgYyAxNDEuMTc2NDA1LDEyLjgxMDc1MDUgMTQxLjA2MTM1OCwxMy4wNDc1Mzk2IDE0MC42NTAwMjIsMTMuNDU4ODc1OSBjIDE0MC42NDI0NTQsMTMuNDY2NDQzOSAxNDAuNDc1MzI2LDEzLjg2NzE0NDUgMTQwLjQ3OTQ0MSwxNC4yNzEwOTI0IGMgMTQwLjQ4MzIzMiwxNC42NDMzMzk1IDE0MC40OTYyNDQsMTUuMDA1NTgzMiAxNDAuNDA5NDI0LDE1LjE3OTIyMzUgYyAxNDAuMTE1NDAxLDE1Ljc2NzI2ODYgMTQwLjA3OTI5MiwxNS41Mjg4OTk3IDEzOS44NDgwMDcsMTUuOTkzNTM0NCBjIDEzOS43NTc5NDEsMTYuMTc0NDcwMSAxMzkuNTgwMzYsMTYuNjg4OTI3OCAxMzkuNjk5MzI0LDE2Ljk4ODk5OTIgYyAxMzkuODg1ODU0LDE3LjQ1OTQ5NzkgMTQwLjM3NDQ3NiwxNy45MDI4ODg2IDE0MC40MDk0MjQsMTguMDc3NjI3NCBjIDE0MC40NzY2OCwxOC40MTM5MDggMTQwLjQ3OTMzNCwxOC42NjI4MzMgMTQwLjQ3OTQ0MSwxOS4xMTcwMzk4IGMgMTQwLjQ3OTU2OSwxOS42NjY1OTc5IDE0MC43MzkwMTksMTkuOTYyMTM4MyAxNDEuMDUzMzcsMjAuMjI2NjI0NSBjIDE0MS4zNDQ2ODYsMjAuNDcxNzI4OCAxNDEuNjgzMTUxLDIwLjY5MDE2MzUgMTQxLjkwNjA3MywyMS4wNTkzNzkyIGMgMTQyLjM2OTU0NywyMS44MjcwMDY3IDE0Mi4yMDU4ODQsMjIuMzMxMTk4NiAxNDIuMzY5NTQ3LDIyLjM5NTM3MzQgYyAxNDIuNzYwNjAzLDIyLjU0ODcxMzMgMTQzLjMyODg2NiwyMi40OTY2NTU2IDE0My42MDA2NDksMjIuNzY4NDM4IGMgMTQzLjg4MzgyMSwyMy4wNTE2MSAxNDMuODMyMzg1LDIzLjk5OTU0MDMgMTQ0LjA5MzA4OSwyNC4wNzI5NDE4IGMgMTQ0Ljc4NjczLDI0LjI2ODIzNjkgMTQ0LjkxODQ1LDI0LjMzMjY2MDQgMTQ1LjI2NjQyNiwyNC40NDg1Mjk1IGMgMTQ1LjYxNDQwMywyNC41NjQzOTg3IDE0NS42ODgwMjUsMjQuOTI2NDg3OSAxNDYuMTY1OTgyLDI1LjMzMjAyNjIgYyAxNDYuNjQzOTM4LDI1LjczNzU2NDYgMTQ2Ljg4ODU4NCwyNi4wMDY0NTMxIDE0Ny4yMDg2NSwyNi4zNzQ4NDEyIGMgMTQ3LjQ3MzU3NiwyNi42Nzk3NjM0IDE0Ny44Mjg5NiwyNi43MDIyNjI2IDE0OC4xNzg1MTMsMjYuNzkyNjUwNyBjIDE0OC40NTMyODMsMjYuODYzNzAxIDE0OC43MjQ0NDksMjYuOTc2Njk4OSAxNDguOTQ1MjQ3LDI3LjMwMTc4ODggYyAxNDkuMTA1MjgxLDI3LjUzNzQxMzcgMTQ5LjU5MjI0LDI3LjU4OTc2MjIgMTUwLjE3NDcwOCwyNy41OTEwNDExIGMgMTUwLjYyNTY4OSwyNy41OTIwMzEzIDE1MS4xMzM5MjcsMjcuNTYyNDA2MiAxNTEuNTkyMDA3LDI3LjU2MzUzMDIgYyAxNTIuMDUxNjgsMjcuNTY0NjU4MiAxNTIuNDYwODQ3LDI3LjU5Njc0OTQgMTUyLjcxMDk3LDI3LjcyMTgxMDggYyAxNTMuMTE5Nzc2LDI3LjkyNjIxMzggMTUyLjk0Nzg2LDI3LjkwNDUyNDUgMTUzLjE5MDE2OSwyOC4xNDY4MzM0IGMgMTUzLjUwMzY4MiwyOC40NjAzNDYgMTU0LjA3MzY1OSwyOC40NTY5ODM5IDE1NC42NTQwNjIsMjguNDM1NzM1MyBjIDE1NS4xODU0NTYsMjguNDE2MjgxIDE1NS43MjU1OSwyOC4zODE4MzMzIDE1Ni4wODU2MzYsMjguNTYxODU2NiBjIDE1Ni40Mzc5NjYsMjguNzM4MDIxMSAxNTYuMzY5ODA1LDI4LjkyODQzNjMgMTU2Ljk0MDE2NiwyOS4yMTM2MTY4IGMgMTU3LjA3NDMyNCwyOS4yODA2OTYgMTU3LjQ5NzI2NCwyOS4yNjc2OTcgMTU3LjkwODg0OCwyOS4yOTI1MTgzIGMgMTU4LjE4NjgwOCwyOS4zMDkyODEyIDE1OC40NTk1ODgsMjkuMzQzMjkzNSAxNTguNjM0NzQxLDI5LjQzMDg2OTYgYyAxNTkuMDY5MjQ4LDI5LjY0ODEyMjQgMTU5LjAzMjU1LDI5Ljg1OTg4OTcgMTU5LjE5NTc1NiwyOS45MTQyOTE5IGMgMTU5LjU3OTc5NywzMC4wNDIzMDU0IDE2MC4zODY2MTYsMzAuMTY2OTY1NCAxNjAuNjA5NTU0LDMwLjM4OTkwMzIgYyAxNjAuOTk4OTYsMzAuNzc5MzA5NSAxNjEuMzU5MTg4LDMxLjI1MDg2MjggMTYxLjY0NzMxOSwzMS41MTY1MDEzIGMgMTYxLjk4NzU4MiwzMS44MzAyMDMxIDE2Mi4yMTIzOCwzMi4wNjE3NTkgMTYyLjMzMjUzNywzMi41MTM1MDMgYyAxNjIuNDAzNzAzLDMyLjc4MTA2MjMgMTYyLjQzODE2MiwzMy4xMjU4NjI5IDE2Mi40MzgxNjIsMzMuNjEwNzIwNyBjIDE2Mi40MzgxNjIsMzQuMDE4NjYxOSAxNjIuNDUyODM3LDM0LjczNzQwODUgMTYyLjIyNjY2MSwzNC45NjM1ODQ1IGMgMTYyLjA5OTQxNSwzNS4wOTA4Mjk2IDE2MS4zNDMwNCwzNS4xNzQwODQ4IDE2MS4xODQ2OTksMzUuMjUzMjU1NSBjIDE2MC44MDMxNTMsMzUuNDQ0MDI4NyAxNjAuNDM2NjM2LDM1Ljg5OTA4MTUgMTYwLjEyNjU0NywzNi4yMDkxNzA0IGMgMTU5Ljg0MTc5NCwzNi40OTM5MjMzIDE1OS4xNzkxNzcsMzcuMjQ0MjI1MSAxNTguODQ5MzQ0LDM3LjUyNzE3MjcgYyAxNTguNTYyMzIzLDM3Ljc3MzM5NDcgMTU3Ljk0NTgzLDM3LjY4ODc4MzEgMTU3LjY0OTg1OSwzNy44NDU4MTI5IGMgMTU3LjI4ODU1OSwzOC4wMzc1MDM1IDE1Ny4zMDQzMzQsMzguMjI5NDIzMiAxNTcuMDI4NzM5LDM4LjM2NzIyMDMgYyAxNTYuNzc5NjE1LDM4LjQ5MTc4MjcgMTU2LjM0OTIzMiwzOC41MzI1NDQ1IDE1NS44Njk5MTUsMzguNTQ0NDY4MSBjIDE1NS40ODI0NDEsMzguNTU0MTA3IDE1NS4wNjI5ODgsMzguNTQ0OTAwNSAxNTQuNjgxNDYsMzguNTQ1ODgzOCBjIDE1NC4xMzE4OSwzOC41NDczMDAyIDE1My42NjEwMSwzOC41Njk4NTkzIDE1My40Nzc3NDIsMzguNzAwMzQwNSBjIDE1My4wODkxNzIsMzguOTc2OTkwNCAxNTMuMDE5MjAyLDM5LjI0NjI5MiAxNTIuNDQ4NzI1LDM5LjMzNTA0MjQgYyAxNTIuMTgxNTM5LDM5LjM3NjYwOTEgMTUxLjgwNDU2NCwzOS4zNzg1NzA4IDE1MS4yMzM2NDgsMzkuMzIzMTMzNCBjIDE1MC40MTg3NDEsMzkuMjQ0MDAzOSAxNTAuNTk2MTA2LDM4Ljg1OTY2MDUgMTUwLjE5Nzk0MSwzOC43MDAzNDA1IGMgMTUwLjAyOTAwNSwzOC42MzI3NDMgMTQ5LjY4MDM4NCwzOC41OTUzNzQxIDE0OS4zMjEwMzMsMzguNTYwODE5OCBjIDE0OC44MzM0MzMsMzguNTEzOTMzNCAxNDguMzI2MDc5LDM4LjQ3MjIyOTIgMTQ4LjIyMTA3LDM4LjM2NzIyMDMgYyAxNDcuOTYyODE0LDM4LjEwODk2NDYgMTQ4LjA0NTg5OCwzOC4wNzY3Mjk5IDE0Ny43ODY1NjMsMzcuOTU0ODM3NSBjIDE0Ny40NzU1MjcsMzcuODA4NjQ0MSAxNDYuOTc3OTQ1LDM3LjcxMjU1MTcgMTQ2LjYwNjcyOSwzNy40Njc3MDc0IGMgMTQ2LjIxODc5NiwzNy4yMTE4MzcyIDE0NS44NDcxNCwzNi44MTA2NDgyIDE0NS40NjUxNjcsMzYuNDAyOTYzMiBjIDE0NS4xNzY4NzEsMzYuMDk1MjYxMyAxNDQuODgyNjk5LDM1Ljc4Mzg1ODkgMTQ0LjU3MTIxNiwzNS41Mjg0NDI5IGMgMTQ0LjMzMDM4NywzNS4zMzA5NjI3IDE0My43MzExNywzNS4xNjkwOTI5IDE0My4zMDQzMzUsMzUuMDM2MDAyNyBjIDE0Mi44Nzc1LDM0LjkwMjkxMjUgMTQzLjA4MDksMzQuMDU3Nzc4NCAxNDIuODczOTM2LDMzLjgxNTQyNzEgYyAxNDIuNjAwNjY3LDMzLjQ5NTQzMzggMTQyLjE0ODM1OSwzMy4zNzkwMjMzIDE0MS40NDY5MTksMzMuNDEzODgxNSBjIDE0MS4wNDc3NCwzMy40MzM3MTg4IDE0MC41Njc4NzgsMzMuNTAyNTQ1OCAxMzkuOTk0NDE1LDMzLjYxMDcyMDcgYyAxMzkuNjI0ODcsMzMuNjgwNDI5NiAxMzkuNjE3ODQzLDM0LjA5NDU3MTQgMTM4Ljk1MTU5OSwzNC41ODg0NjkxIGwgMTM4LjgxMTUxNywzNS4xMDYxMDExIHogXCI7XG4vLyB2YXIgcGF0aCA9IHBvaW50MSArIHBvaW50MiArIHBvaW50MyArIHBvaW50NCArIHBvaW50NSArIHBvaW50NiArIHBvaW50NyArIHBvaW50OCArIHBvaW50OTtcbnZhciBwYXRoID0gbGV0dHJlc0FOVCArIGxldHRyZVM7XG5cbm1vZHVsZS5leHBvcnRzID0gcGF0aDsiLCIndXNlIHN0cmljdCc7XG5cbnZhciBzcXJ0ID0gTWF0aC5zcXJ0O1xudmFyIHBvdyA9IE1hdGgucG93O1xuXG5mdW5jdGlvbiBzaWduKHgpIHtcblx0cmV0dXJuIHggPyB4IDwgMCA/IC0xIDogMSA6IDA7XG59XG5cbmZ1bmN0aW9uIHJhbmdlKHN0YXJ0LCBjb3VudCkge1xuICAgIHJldHVybiBBcnJheS5hcHBseSgwLCBBcnJheShjb3VudCkpLm1hcChmdW5jdGlvbiAoZWxlbWVudCwgaW5kZXgpIHtcbiAgICBcdHJldHVybiBpbmRleCArIHN0YXJ0XG4gICAgfSk7XG59XG5cbmZ1bmN0aW9uIGRpc3RhbmNlKGEsIGIpe1xuXHRyZXR1cm4gc3FydChwb3coYS54IC0gYi54LCAyKSArIHBvdyhhLnkgLSBiLnksIDIpKTtcbn1cblxuZnVuY3Rpb24gbm9ybSh2KXtcblx0cmV0dXJuIHNxcnQocG93KHYueCwgMikgKyBwb3codi55LCAyKSk7XG59XG5cbm1vZHVsZS5leHBvcnRzID0ge1xuXHRzaWduOiBzaWduLFxuXHRyYW5nZTogcmFuZ2UsXG5cdGRpc3RhbmNlOiBkaXN0YW5jZSxcblx0bm9ybTogbm9ybVxufSIsIid1c2Ugc3RyaWN0J1xuXG5mdW5jdGlvbiBWZWN0b3IoeCwgeSkge1xuICAgIHRoaXMueCA9IHg7ICAgICAgICAgICAgICAgIFxuICAgIHRoaXMueSA9IHk7XG59XG5cblZlY3Rvci5wcm90b3R5cGUubm9ybSA9IGZ1bmN0aW9uKCl7XG5cdHJldHVybiBNYXRoLnNxcnQodGhpcy54ICogdGhpcy54ICsgdGhpcy55ICogdGhpcy55KTtcbn1cblxuVmVjdG9yLnByb3RvdHlwZS5ub3JtYWxpemUgPSBmdW5jdGlvbigpe1xuXHR2YXIgbm9ybSA9IHRoaXMubm9ybSgpO1xuXHR0aGlzLnggPSB0aGlzLnggLyBub3JtO1xuXHR0aGlzLnkgPSB0aGlzLnkgLyBub3JtO1xufVxuXG5cblxubW9kdWxlLmV4cG9ydHMgPSBWZWN0b3I7IiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgX2FudENvbG9ueSA9IHJlcXVpcmUoJy4vaW5kZXguanMnKTtcblxudmFyIGNvbnRhaW5lciA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3IoJy5jb2xvbnknKTtcblxudmFyIG9wdGlvbnMgPSB7XG5cdHZlbG9jaXR5OiAwLjAwMSxcblx0bmJBbnRzOiA0MDAwLFxuXHR3ZWlnaHQ6IDEwLFxuXHRyZXBTaXplOiAwLjA1LFxuXHRyZXBTcGVlZDogMC4wMDIsXG5cdG5iU3RhcnQ6IDUwMCxcblx0bmJSYW5kOiA5MDBcblx0Ly8gb2JqIHBhciBkZWZhdXRcbn07XG5cbnZhciBhbnRDb2xvbnkgPSBfYW50Q29sb255KGNvbnRhaW5lciwgb3B0aW9ucyk7XG5cbndpbmRvdy5hZGRFdmVudExpc3RlbmVyKCdjbGljaycsIGZ1bmN0aW9uICgpe1xuXHQvLyBvcHRpb25zLnZlbG9jaXR5ID0gMC4wMDM7XG5cdG9wdGlvbnMubmJBbnRzID0gMjAwMDA7XG5cdC8vIG9wdGlvbnMud2VpZ2h0ID0gMTAwMDAwMDA7XG5cdC8vIG9wdGlvbnMucmVwU3BlZWQgPSAwLjAxO1xuXHQvLyBvcHRpb25zLnJlcFNpemUgPSAwLjE7XG5cblx0Ly8gYW50Q29sb255LmNoYW5nZU9wdGlvbnMob3B0aW9ucyk7XG5cdGFudENvbG9ueS5jaGFuZ2VPcHRpb25zKG9wdGlvbnMpO1xufSk7XG5cbiJdfQ==
