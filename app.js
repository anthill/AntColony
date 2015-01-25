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
var WEIGHT = 10;
// var REPULSION = 0.05;
// var REPULSIONSPEED = 0.002;
// var ANTVELOCITY = 0.001;

module.exports = function(container, initPoints, options){

    console.log('Options ant :', options);
    // Define those parameters as attributes of Ant object ?
    var REPULSION = options.repSize;
    var REPULSIONSPEED = options.repSpeed;
    var ANTVELOCITY = options.velocity;
    var INTELLIGENCE = options.intelligence;

    var mouse = liveMousePosition(container);

    var points = initPoints.points;
    var citySet = initPoints.citySet;
    var textPointsId = initPoints.textPointsId;
    var possibleStartPointsId = initPoints.possibleStartPointsId;


    function Ant(point) {
        this.x = point.x;                
        this.y = point.y;
        this.velocity = ANTVELOCITY;
        this.intelligence = INTELLIGENCE;
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
                        this.edges.forEach(function(e){
                            var a = e.pt1, b = e.pt2, weight = 1;  
                            // increased dropped pheromons for textEdges
                            if ((citySet.indexOf(a.id) != -1) && citySet.indexOf(b.id) != -1 && (Math.abs(a.id - b.id) == 1))
                            {
                                weight *= WEIGHT;
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
            if (random() < this.intelligence){
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
	var objPopulationMax = options.nbAnts;
	var objPopulationAdjusted = objPopulationMax;
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
		if (antNumber < objPopulationAdjusted - 50){
			FPSMonitor.style.color = "green";
			population = antsGroup.create(population);
		}	
		else if (antNumber > objPopulationAdjusted){
			population = antsGroup.remove(population, antNumber - objPopulationAdjusted);
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
			objPopulationAdjusted = objPopulationAdjusted * maxDeltaTime / deltaTime;
			FPSOverLimitCount = 0;
		}

		while (FPSUnderLimitCount > 50 && objPopulationAdjusted < objPopulationMax) {
			objPopulationAdjusted += 10;
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
		objPopulationMax = opts.nbAnts;
		objPopulationAdjusted = objPopulationMax;

		population.forEach(function(ant){
			ant.velocity = opts.velocity;
			ant.intelligence = opts.intelligence;
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

},{"./ant.js":"/Users/Romain/Documents/Programmation/ants/AntColony/src/ant.js","./antsGroup":"/Users/Romain/Documents/Programmation/ants/AntColony/src/antsGroup.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/runAnimation.js":[function(require,module,exports){
'use strict';

require('../index.js')(document.body.querySelector('.colony'));
},{"../index.js":"/Users/Romain/Documents/Programmation/ants/AntColony/index.js"}],"/Users/Romain/Documents/Programmation/ants/AntColony/src/svgPath.js":[function(require,module,exports){
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
},{}]},{},["/Users/Romain/Documents/Programmation/ants/AntColony/src/runAnimation.js"])
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi91c3IvbG9jYWwvbGliL25vZGVfbW9kdWxlcy93YXRjaGlmeS9ub2RlX21vZHVsZXMvYnJvd3NlcmlmeS9ub2RlX21vZHVsZXMvYnJvd3Nlci1wYWNrL19wcmVsdWRlLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9pbmRleC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9pY2guanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL25vZGVfbW9kdWxlcy90d28tc3VtL3R3by1zdW0uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL3JvYnVzdC1zY2FsZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL25vZGVfbW9kdWxlcy9yb2J1c3Qtc3VidHJhY3Qvcm9idXN0LWRpZmYuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1bS9yb2J1c3Qtc3VtLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vbm9kZV9tb2R1bGVzL3R3by1wcm9kdWN0L3R3by1wcm9kdWN0LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9yb2J1c3Qtb3JpZW50YXRpb24vb3JpZW50YXRpb24uanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3NpbXBsaWNpYWwtY29tcGxleC9ub2RlX21vZHVsZXMvYml0LXR3aWRkbGUvdHdpZGRsZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvc2ltcGxpY2lhbC1jb21wbGV4L25vZGVfbW9kdWxlcy91bmlvbi1maW5kL2luZGV4LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvdG9wb2xvZ3kuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvdW5pcS91bmlxLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvdHJpYW5ndWxhdGUuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L25vZGVfbW9kdWxlcy9wYXJzZS1zdmctcGF0aC9pbmRleC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2FudC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2FudHNHcm91cC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2NyZWF0ZUVkZ2VzLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvZWRnZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL2luaXRpYWxpemVQb2ludHMuanMiLCIvVXNlcnMvUm9tYWluL0RvY3VtZW50cy9Qcm9ncmFtbWF0aW9uL2FudHMvQW50Q29sb255L3NyYy9tb3VzZS5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3BvaW50LmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvcmVuZGVyaW5nLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvcnVuQW5pbWF0aW9uLmpzIiwiL1VzZXJzL1JvbWFpbi9Eb2N1bWVudHMvUHJvZ3JhbW1hdGlvbi9hbnRzL0FudENvbG9ueS9zcmMvc3ZnUGF0aC5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3V0aWxpdGllcy5qcyIsIi9Vc2Vycy9Sb21haW4vRG9jdW1lbnRzL1Byb2dyYW1tYXRpb24vYW50cy9BbnRDb2xvbnkvc3JjL3ZlY3Rvci5qcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtBQ0FBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3ZDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzdiQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2pEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDM0pBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMzSkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM3TEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNU1BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3pEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN0VkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDekRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM5SkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN2REE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNuTkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ25DQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3JGQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM5RUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDeEtBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDbkJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ1RBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdlBBO0FBQ0E7QUFDQTs7QUNGQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzVCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gZSh0LG4scil7ZnVuY3Rpb24gcyhvLHUpe2lmKCFuW29dKXtpZighdFtvXSl7dmFyIGE9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtpZighdSYmYSlyZXR1cm4gYShvLCEwKTtpZihpKXJldHVybiBpKG8sITApO3ZhciBmPW5ldyBFcnJvcihcIkNhbm5vdCBmaW5kIG1vZHVsZSAnXCIrbytcIidcIik7dGhyb3cgZi5jb2RlPVwiTU9EVUxFX05PVF9GT1VORFwiLGZ9dmFyIGw9bltvXT17ZXhwb3J0czp7fX07dFtvXVswXS5jYWxsKGwuZXhwb3J0cyxmdW5jdGlvbihlKXt2YXIgbj10W29dWzFdW2VdO3JldHVybiBzKG4/bjplKX0sbCxsLmV4cG9ydHMsZSx0LG4scil9cmV0dXJuIG5bb10uZXhwb3J0c312YXIgaT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2Zvcih2YXIgbz0wO288ci5sZW5ndGg7bysrKXMocltvXSk7cmV0dXJuIHN9KSIsIid1c2Ugc3RyaWN0JztcblxudmFyIGluaXRSZW5kZXJpbmcgPSByZXF1aXJlKCcuL3NyYy9yZW5kZXJpbmcuanMnKTtcbnZhciBpbml0aWFsaXplUG9pbnRzID0gcmVxdWlyZSgnLi9zcmMvaW5pdGlhbGl6ZVBvaW50cy5qcycpO1xudmFyIGNyZWF0ZUVkZ2VzID0gcmVxdWlyZSgnLi9zcmMvY3JlYXRlRWRnZXMuanMnKTtcbi8vIHZhciBpbml0QW50cyA9IHJlcXVpcmUoJy4vc3JjL2luaXRpYWxpemVBbnRzJyk7XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24gaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRpb25zKXtcblxuXHR2YXIgcmVuZGVyLCBwb2ludHNJbmZvcywgZWRnZXMsIHBvcHVsYXRpb24sIHBvaW50c01hcDtcblxuXG5cdGZ1bmN0aW9uIF9pbml0KGNvbnRhaW5lckVsZW1lbnQsIG9wdGlvbnMpe1xuXHRcdHBvaW50c0luZm9zID0gaW5pdGlhbGl6ZVBvaW50cyhvcHRpb25zLm5iU3RhcnQsIG9wdGlvbnMubmJSYW5kKTtcblx0XHRlZGdlcyA9IGNyZWF0ZUVkZ2VzKHBvaW50c0luZm9zLnBvaW50cyk7XG5cdFx0Ly8gcG9wdWxhdGlvbiA9IG9wdGlvbnMubmJBbnRzO1xuXHRcdC8vIHBvcHVsYXRpb24gPSBpbml0QW50cyhjb250YWluZXJFbGVtZW50LCBwb2ludHNJbmZvcywgb3B0aW9ucyk7XG5cdFx0cG9pbnRzTWFwID0ge1xuXHRcdFx0cG9pbnRzSW5mb3M6IHBvaW50c0luZm9zLFxuXHRcdFx0ZWRnZXM6IGVkZ2VzXG5cdFx0XHQvLyBwb3B1bGF0aW9uOiBwb3B1bGF0aW9uXG5cdFx0fTtcblx0XHRyZW5kZXIgPSBpbml0UmVuZGVyaW5nKGNvbnRhaW5lckVsZW1lbnQsIHBvaW50c01hcCwgb3B0aW9ucyk7XG5cdH1cblxuXHRfaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRpb25zKTtcblxuXHRyZXR1cm4ge1xuXHRcdHRvZ2dsZVBsYXlQYXVzZTogZnVuY3Rpb24oKXsgcmVuZGVyLnRvZ2dsZVBsYXlQYXVzZSgpIH0sXG5cdFx0Y2hhbmdlT3B0aW9uczogZnVuY3Rpb24ob3B0cyl7XG5cdFx0XHRyZW5kZXIubW9kaWZ5QW50cyhvcHRzKTtcblx0XHR9LFxuXHRcdHJlc2V0OiBmdW5jdGlvbihvcHRzKXtcblx0XHRcdHJlbmRlci5yZXNldCgpO1xuXG5cdFx0XHRcdC8vIHJlc2V0IGVsZW1lbnRzXG5cdFx0XHRfaW5pdChjb250YWluZXJFbGVtZW50LCBvcHRzKTtcblx0XHR9XG5cdH07XG59OyIsIlwidXNlIHN0cmljdFwiXG5cbi8vSGlnaCBsZXZlbCBpZGVhOlxuLy8gMS4gVXNlIENsYXJrc29uJ3MgaW5jcmVtZW50YWwgY29uc3RydWN0aW9uIHRvIGZpbmQgY29udmV4IGh1bGxcbi8vIDIuIFBvaW50IGxvY2F0aW9uIGluIHRyaWFuZ3VsYXRpb24gYnkganVtcCBhbmQgd2Fsa1xuXG5tb2R1bGUuZXhwb3J0cyA9IGluY3JlbWVudGFsQ29udmV4SHVsbFxuXG52YXIgb3JpZW50ID0gcmVxdWlyZShcInJvYnVzdC1vcmllbnRhdGlvblwiKVxudmFyIGNvbXBhcmVDZWxsID0gcmVxdWlyZShcInNpbXBsaWNpYWwtY29tcGxleFwiKS5jb21wYXJlQ2VsbHNcblxuZnVuY3Rpb24gY29tcGFyZUludChhLCBiKSB7XG4gIHJldHVybiBhIC0gYlxufVxuXG5mdW5jdGlvbiBTaW1wbGV4KHZlcnRpY2VzLCBhZGphY2VudCwgYm91bmRhcnkpIHtcbiAgdGhpcy52ZXJ0aWNlcyA9IHZlcnRpY2VzXG4gIHRoaXMuYWRqYWNlbnQgPSBhZGphY2VudFxuICB0aGlzLmJvdW5kYXJ5ID0gYm91bmRhcnlcbiAgdGhpcy5sYXN0VmlzaXRlZCA9IC0xXG59XG5cblNpbXBsZXgucHJvdG90eXBlLmZsaXAgPSBmdW5jdGlvbigpIHtcbiAgdmFyIHQgPSB0aGlzLnZlcnRpY2VzWzBdXG4gIHRoaXMudmVydGljZXNbMF0gPSB0aGlzLnZlcnRpY2VzWzFdXG4gIHRoaXMudmVydGljZXNbMV0gPSB0XG4gIHZhciB1ID0gdGhpcy5hZGphY2VudFswXVxuICB0aGlzLmFkamFjZW50WzBdID0gdGhpcy5hZGphY2VudFsxXVxuICB0aGlzLmFkamFjZW50WzFdID0gdVxufVxuXG5mdW5jdGlvbiBHbHVlRmFjZXQodmVydGljZXMsIGNlbGwsIGluZGV4KSB7XG4gIHRoaXMudmVydGljZXMgPSB2ZXJ0aWNlc1xuICB0aGlzLmNlbGwgPSBjZWxsXG4gIHRoaXMuaW5kZXggPSBpbmRleFxufVxuXG5mdW5jdGlvbiBjb21wYXJlR2x1ZShhLCBiKSB7XG4gIHJldHVybiBjb21wYXJlQ2VsbChhLnZlcnRpY2VzLCBiLnZlcnRpY2VzKVxufVxuXG5mdW5jdGlvbiBiYWtlT3JpZW50KGQpIHtcbiAgdmFyIGNvZGUgPSBbXCJmdW5jdGlvbiBvcmllbnQoKXt2YXIgdHVwbGU9dGhpcy50dXBsZTtyZXR1cm4gdGVzdChcIl1cbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIGlmKGkgPiAwKSB7XG4gICAgICBjb2RlLnB1c2goXCIsXCIpXG4gICAgfVxuICAgIGNvZGUucHVzaChcInR1cGxlW1wiLCBpLCBcIl1cIilcbiAgfVxuICBjb2RlLnB1c2goXCIpfXJldHVybiBvcmllbnRcIilcbiAgdmFyIHByb2MgPSBuZXcgRnVuY3Rpb24oXCJ0ZXN0XCIsIGNvZGUuam9pbihcIlwiKSlcbiAgdmFyIHRlc3QgPSBvcmllbnRbZCsxXVxuICBpZighdGVzdCkge1xuICAgIHRlc3QgPSBvcmllbnRcbiAgfVxuICByZXR1cm4gcHJvYyh0ZXN0KVxufVxuXG52YXIgQkFLRUQgPSBbXVxuXG5mdW5jdGlvbiBUcmlhbmd1bGF0aW9uKGRpbWVuc2lvbiwgdmVydGljZXMsIHNpbXBsaWNlcykge1xuICB0aGlzLmRpbWVuc2lvbiA9IGRpbWVuc2lvblxuICB0aGlzLnZlcnRpY2VzID0gdmVydGljZXNcbiAgdGhpcy5zaW1wbGljZXMgPSBzaW1wbGljZXNcbiAgdGhpcy5pbnRlcmlvciA9IHNpbXBsaWNlcy5maWx0ZXIoZnVuY3Rpb24oYykge1xuICAgIHJldHVybiAhYy5ib3VuZGFyeVxuICB9KVxuXG4gIHRoaXMudHVwbGUgPSBuZXcgQXJyYXkoZGltZW5zaW9uKzEpXG4gIGZvcih2YXIgaT0wOyBpPD1kaW1lbnNpb247ICsraSkge1xuICAgIHRoaXMudHVwbGVbaV0gPSB0aGlzLnZlcnRpY2VzW2ldXG4gIH1cblxuICB2YXIgbyA9IEJBS0VEW2RpbWVuc2lvbl1cbiAgaWYoIW8pIHtcbiAgICBvID0gQkFLRURbZGltZW5zaW9uXSA9IGJha2VPcmllbnQoZGltZW5zaW9uKVxuICB9XG4gIHRoaXMub3JpZW50ID0gb1xufVxuXG52YXIgcHJvdG8gPSBUcmlhbmd1bGF0aW9uLnByb3RvdHlwZVxuXG4vL0RlZ2VuZXJhdGUgc2l0dWF0aW9uIHdoZXJlIHdlIGFyZSBvbiBib3VuZGFyeSwgYnV0IGNvcGxhbmFyIHRvIGZhY2VcbnByb3RvLmhhbmRsZUJvdW5kYXJ5RGVnZW5lcmFjeSA9IGZ1bmN0aW9uKGNlbGwsIHBvaW50KSB7XG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIG4gPSB0aGlzLnZlcnRpY2VzLmxlbmd0aCAtIDFcbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG5cbiAgLy9EdW1iIHNvbHV0aW9uOiBKdXN0IGRvIGRmcyBmcm9tIGJvdW5kYXJ5IGNlbGwgdW50aWwgd2UgZmluZCBhbnkgcGVhaywgb3IgdGVybWluYXRlXG4gIHZhciB0b1Zpc2l0ID0gWyBjZWxsIF1cbiAgY2VsbC5sYXN0VmlzaXRlZCA9IC1uXG4gIHdoaWxlKHRvVmlzaXQubGVuZ3RoID4gMCkge1xuICAgIGNlbGwgPSB0b1Zpc2l0LnBvcCgpXG4gICAgdmFyIGNlbGxWZXJ0cyA9IGNlbGwudmVydGljZXNcbiAgICB2YXIgY2VsbEFkaiA9IGNlbGwuYWRqYWNlbnRcbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICB2YXIgbmVpZ2hib3IgPSBjZWxsQWRqW2ldXG4gICAgICBpZighbmVpZ2hib3IuYm91bmRhcnkgfHwgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPD0gLW4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHZhciBudiA9IG5laWdoYm9yLnZlcnRpY2VzXG4gICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgIHZhciB2diA9IG52W2pdXG4gICAgICAgIGlmKHZ2IDwgMCkge1xuICAgICAgICAgIHR1cGxlW2pdID0gcG9pbnRcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICB0dXBsZVtqXSA9IHZlcnRzW3Z2XVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICB2YXIgbyA9IHRoaXMub3JpZW50KClcbiAgICAgIGlmKG8gPiAwKSB7XG4gICAgICAgIHJldHVybiBuZWlnaGJvclxuICAgICAgfVxuICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSAtblxuICAgICAgaWYobyA9PT0gMCkge1xuICAgICAgICB0b1Zpc2l0LnB1c2gobmVpZ2hib3IpXG4gICAgICB9XG4gICAgfVxuICB9XG4gIHJldHVybiBudWxsXG59XG5cbnByb3RvLndhbGsgPSBmdW5jdGlvbihwb2ludCwgcmFuZG9tKSB7XG4gIC8vQWxpYXMgbG9jYWwgcHJvcGVydGllc1xuICB2YXIgbiA9IHRoaXMudmVydGljZXMubGVuZ3RoIC0gMVxuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuXG4gIC8vQ29tcHV0ZSBpbml0aWFsIGp1bXAgY2VsbFxuICB2YXIgaW5pdEluZGV4ID0gcmFuZG9tID8gKHRoaXMuaW50ZXJpb3IubGVuZ3RoICogTWF0aC5yYW5kb20oKSl8MCA6ICh0aGlzLmludGVyaW9yLmxlbmd0aC0xKVxuICB2YXIgY2VsbCA9IHRoaXMuaW50ZXJpb3JbIGluaXRJbmRleCBdXG5cbiAgLy9TdGFydCB3YWxraW5nXG5vdXRlckxvb3A6XG4gIHdoaWxlKCFjZWxsLmJvdW5kYXJ5KSB7XG4gICAgdmFyIGNlbGxWZXJ0cyA9IGNlbGwudmVydGljZXNcbiAgICB2YXIgY2VsbEFkaiA9IGNlbGwuYWRqYWNlbnRcblxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIHR1cGxlW2ldID0gdmVydHNbY2VsbFZlcnRzW2ldXVxuICAgIH1cbiAgICBjZWxsLmxhc3RWaXNpdGVkID0gblxuXG4gICAgLy9GaW5kIGZhcnRoZXN0IGFkamFjZW50IGNlbGxcbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICB2YXIgbmVpZ2hib3IgPSBjZWxsQWRqW2ldXG4gICAgICBpZihuZWlnaGJvci5sYXN0VmlzaXRlZCA+PSBuKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICB2YXIgcHJldiA9IHR1cGxlW2ldXG4gICAgICB0dXBsZVtpXSA9IHBvaW50XG4gICAgICB2YXIgbyA9IHRoaXMub3JpZW50KClcbiAgICAgIHR1cGxlW2ldID0gcHJldlxuICAgICAgaWYobyA8IDApIHtcbiAgICAgICAgY2VsbCA9IG5laWdoYm9yXG4gICAgICAgIGNvbnRpbnVlIG91dGVyTG9vcFxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgaWYoIW5laWdoYm9yLmJvdW5kYXJ5KSB7XG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSBuXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSAtblxuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICAgIHJldHVyblxuICB9XG5cbiAgcmV0dXJuIGNlbGxcbn1cblxucHJvdG8uYWRkUGVha3MgPSBmdW5jdGlvbihwb2ludCwgY2VsbCkge1xuICB2YXIgbiA9IHRoaXMudmVydGljZXMubGVuZ3RoIC0gMVxuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcbiAgdmFyIHR1cGxlID0gdGhpcy50dXBsZVxuICB2YXIgaW50ZXJpb3IgPSB0aGlzLmludGVyaW9yXG4gIHZhciBzaW1wbGljZXMgPSB0aGlzLnNpbXBsaWNlc1xuXG4gIC8vV2Fsa2luZyBmaW5pc2hlZCBhdCBib3VuZGFyeSwgdGltZSB0byBhZGQgcGVha3NcbiAgdmFyIHRvdmlzaXQgPSBbIGNlbGwgXVxuXG4gIC8vU3RyZXRjaCBpbml0aWFsIGJvdW5kYXJ5IGNlbGwgaW50byBhIHBlYWtcbiAgY2VsbC5sYXN0VmlzaXRlZCA9IG5cbiAgY2VsbC52ZXJ0aWNlc1tjZWxsLnZlcnRpY2VzLmluZGV4T2YoLTEpXSA9IG5cbiAgY2VsbC5ib3VuZGFyeSA9IGZhbHNlXG4gIGludGVyaW9yLnB1c2goY2VsbClcblxuICAvL1JlY29yZCBhIGxpc3Qgb2YgYWxsIG5ldyBib3VuZGFyaWVzIGNyZWF0ZWQgYnkgYWRkZWQgcGVha3Mgc28gd2UgY2FuIGdsdWUgdGhlbSB0b2dldGhlciB3aGVuIHdlIGFyZSBhbGwgZG9uZVxuICB2YXIgZ2x1ZUZhY2V0cyA9IFtdXG5cbiAgLy9EbyBhIHRyYXZlcnNhbCBvZiB0aGUgYm91bmRhcnkgd2Fsa2luZyBvdXR3YXJkIGZyb20gc3RhcnRpbmcgcGVha1xuICB3aGlsZSh0b3Zpc2l0Lmxlbmd0aCA+IDApIHtcbiAgICAvL1BvcCBvZmYgcGVhayBhbmQgd2FsayBvdmVyIGFkamFjZW50IGNlbGxzXG4gICAgdmFyIGNlbGwgPSB0b3Zpc2l0LnBvcCgpXG4gICAgdmFyIGNlbGxWZXJ0cyA9IGNlbGwudmVydGljZXNcbiAgICB2YXIgY2VsbEFkaiA9IGNlbGwuYWRqYWNlbnRcbiAgICB2YXIgaW5kZXhPZk4gPSBjZWxsVmVydHMuaW5kZXhPZihuKVxuICAgIGlmKGluZGV4T2ZOIDwgMCkge1xuICAgICAgY29udGludWVcbiAgICB9XG5cbiAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICBpZihpID09PSBpbmRleE9mTikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuXG4gICAgICAvL0ZvciBlYWNoIGJvdW5kYXJ5IG5laWdoYm9yIG9mIHRoZSBjZWxsXG4gICAgICB2YXIgbmVpZ2hib3IgPSBjZWxsQWRqW2ldXG4gICAgICBpZighbmVpZ2hib3IuYm91bmRhcnkgfHwgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPj0gbikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuXG4gICAgICB2YXIgbnYgPSBuZWlnaGJvci52ZXJ0aWNlc1xuXG4gICAgICAvL1Rlc3QgaWYgbmVpZ2hib3IgaXMgYSBwZWFrXG4gICAgICBpZihuZWlnaGJvci5sYXN0VmlzaXRlZCAhPT0gLW4pIHsgICAgICBcbiAgICAgICAgLy9Db21wdXRlIG9yaWVudGF0aW9uIG9mIHAgcmVsYXRpdmUgdG8gZWFjaCBib3VuZGFyeSBwZWFrXG4gICAgICAgIHZhciBpbmRleE9mTmVnMSA9IDBcbiAgICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICAgIGlmKG52W2pdIDwgMCkge1xuICAgICAgICAgICAgaW5kZXhPZk5lZzEgPSBqXG4gICAgICAgICAgICB0dXBsZVtqXSA9IHBvaW50XG4gICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHR1cGxlW2pdID0gdmVydHNbbnZbal1dXG4gICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIHZhciBvID0gdGhpcy5vcmllbnQoKVxuXG4gICAgICAgIC8vVGVzdCBpZiBuZWlnaGJvciBjZWxsIGlzIGFsc28gYSBwZWFrXG4gICAgICAgIGlmKG8gPiAwKSB7XG4gICAgICAgICAgbnZbaW5kZXhPZk5lZzFdID0gblxuICAgICAgICAgIG5laWdoYm9yLmJvdW5kYXJ5ID0gZmFsc2VcbiAgICAgICAgICBpbnRlcmlvci5wdXNoKG5laWdoYm9yKVxuICAgICAgICAgIHRvdmlzaXQucHVzaChuZWlnaGJvcilcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IG5cbiAgICAgICAgICBjb250aW51ZVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gLW5cbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgICB2YXIgbmEgPSBuZWlnaGJvci5hZGphY2VudFxuXG4gICAgICAvL090aGVyd2lzZSwgcmVwbGFjZSBuZWlnaGJvciB3aXRoIG5ldyBmYWNlXG4gICAgICB2YXIgdnZlcnRzID0gY2VsbFZlcnRzLnNsaWNlKClcbiAgICAgIHZhciB2YWRqID0gY2VsbEFkai5zbGljZSgpXG4gICAgICB2YXIgbmNlbGwgPSBuZXcgU2ltcGxleCh2dmVydHMsIHZhZGosIHRydWUpXG4gICAgICBzaW1wbGljZXMucHVzaChuY2VsbClcblxuICAgICAgLy9Db25uZWN0IHRvIG5laWdoYm9yXG4gICAgICB2YXIgb3Bwb3NpdGUgPSBuYS5pbmRleE9mKGNlbGwpXG4gICAgICBpZihvcHBvc2l0ZSA8IDApIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIG5hW29wcG9zaXRlXSA9IG5jZWxsXG4gICAgICB2YWRqW2luZGV4T2ZOXSA9IG5laWdoYm9yXG5cbiAgICAgIC8vQ29ubmVjdCB0byBjZWxsXG4gICAgICB2dmVydHNbaV0gPSAtMVxuICAgICAgdmFkaltpXSA9IGNlbGxcbiAgICAgIGNlbGxBZGpbaV0gPSBuY2VsbFxuXG4gICAgICAvL0ZsaXAgZmFjZXRcbiAgICAgIG5jZWxsLmZsaXAoKVxuXG4gICAgICAvL0FkZCB0byBnbHVlIGxpc3RcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgdmFyIHV1ID0gdnZlcnRzW2pdXG4gICAgICAgIGlmKHV1IDwgMCB8fCB1dSA9PT0gbikge1xuICAgICAgICAgIGNvbnRpbnVlXG4gICAgICAgIH1cbiAgICAgICAgdmFyIG5mYWNlID0gbmV3IEFycmF5KGQtMSlcbiAgICAgICAgdmFyIG5wdHIgPSAwXG4gICAgICAgIGZvcih2YXIgaz0wOyBrPD1kOyArK2spIHtcbiAgICAgICAgICB2YXIgdnYgPSB2dmVydHNba11cbiAgICAgICAgICBpZih2diA8IDAgfHwgayA9PT0gaikge1xuICAgICAgICAgICAgY29udGludWVcbiAgICAgICAgICB9XG4gICAgICAgICAgbmZhY2VbbnB0cisrXSA9IHZ2XG4gICAgICAgIH1cbiAgICAgICAgZ2x1ZUZhY2V0cy5wdXNoKG5ldyBHbHVlRmFjZXQobmZhY2UsIG5jZWxsLCBqKSlcbiAgICAgIH1cbiAgICB9XG4gIH1cblxuICAvL0dsdWUgYm91bmRhcnkgZmFjZXRzIHRvZ2V0aGVyXG4gIGdsdWVGYWNldHMuc29ydChjb21wYXJlR2x1ZSlcblxuICBmb3IodmFyIGk9MDsgaSsxPGdsdWVGYWNldHMubGVuZ3RoOyBpKz0yKSB7XG4gICAgdmFyIGEgPSBnbHVlRmFjZXRzW2ldXG4gICAgdmFyIGIgPSBnbHVlRmFjZXRzW2krMV1cbiAgICB2YXIgYWkgPSBhLmluZGV4XG4gICAgdmFyIGJpID0gYi5pbmRleFxuICAgIGlmKGFpIDwgMCB8fCBiaSA8IDApIHtcbiAgICAgIGNvbnRpbnVlXG4gICAgfVxuICAgIGEuY2VsbC5hZGphY2VudFthLmluZGV4XSA9IGIuY2VsbFxuICAgIGIuY2VsbC5hZGphY2VudFtiLmluZGV4XSA9IGEuY2VsbFxuICB9XG59XG5cbnByb3RvLmluc2VydCA9IGZ1bmN0aW9uKHBvaW50LCByYW5kb20pIHtcbiAgLy9BZGQgcG9pbnRcbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuICB2ZXJ0cy5wdXNoKHBvaW50KVxuXG4gIHZhciBjZWxsID0gdGhpcy53YWxrKHBvaW50LCByYW5kb20pXG4gIGlmKCFjZWxsKSB7XG4gICAgcmV0dXJuXG4gIH1cblxuICAvL0FsaWFzIGxvY2FsIHByb3BlcnRpZXNcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG5cbiAgLy9EZWdlbmVyYXRlIGNhc2U6IElmIHBvaW50IGlzIGNvcGxhbmFyIHRvIGNlbGwsIHRoZW4gd2FsayB1bnRpbCB3ZSBmaW5kIGEgbm9uLWRlZ2VuZXJhdGUgYm91bmRhcnlcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB2diA9IGNlbGwudmVydGljZXNbaV1cbiAgICBpZih2diA8IDApIHtcbiAgICAgIHR1cGxlW2ldID0gcG9pbnRcbiAgICB9IGVsc2Uge1xuICAgICAgdHVwbGVbaV0gPSB2ZXJ0c1t2dl1cbiAgICB9XG4gIH1cbiAgdmFyIG8gPSB0aGlzLm9yaWVudCh0dXBsZSlcbiAgaWYobyA8IDApIHtcbiAgICByZXR1cm5cbiAgfSBlbHNlIGlmKG8gPT09IDApIHtcbiAgICBjZWxsID0gdGhpcy5oYW5kbGVCb3VuZGFyeURlZ2VuZXJhY3koY2VsbCwgcG9pbnQpXG4gICAgaWYoIWNlbGwpIHtcbiAgICAgIHJldHVyblxuICAgIH1cbiAgfVxuXG4gIC8vQWRkIHBlYWtzXG4gIHRoaXMuYWRkUGVha3MocG9pbnQsIGNlbGwpXG59XG5cbi8vRXh0cmFjdCBhbGwgYm91bmRhcnkgY2VsbHNcbnByb3RvLmJvdW5kYXJ5ID0gZnVuY3Rpb24oKSB7XG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIGJvdW5kYXJ5ID0gW11cbiAgdmFyIGNlbGxzID0gdGhpcy5zaW1wbGljZXNcbiAgdmFyIG5jID0gY2VsbHMubGVuZ3RoXG4gIGZvcih2YXIgaT0wOyBpPG5jOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgaWYoYy5ib3VuZGFyeSkge1xuICAgICAgdmFyIGJjZWxsID0gbmV3IEFycmF5KGQpXG4gICAgICB2YXIgY3YgPSBjLnZlcnRpY2VzXG4gICAgICB2YXIgcHRyID0gMFxuICAgICAgdmFyIHBhcml0eSA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgaWYoY3Zbal0gPj0gMCkge1xuICAgICAgICAgIGJjZWxsW3B0cisrXSA9IGN2W2pdXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgcGFyaXR5ID0gaiYxXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIGlmKHBhcml0eSA9PT0gKGQmMSkpIHtcbiAgICAgICAgdmFyIHQgPSBiY2VsbFswXVxuICAgICAgICBiY2VsbFswXSA9IGJjZWxsWzFdXG4gICAgICAgIGJjZWxsWzFdID0gdFxuICAgICAgfVxuICAgICAgYm91bmRhcnkucHVzaChiY2VsbClcbiAgICB9XG4gIH1cbiAgcmV0dXJuIGJvdW5kYXJ5XG59XG5cbmZ1bmN0aW9uIGluY3JlbWVudGFsQ29udmV4SHVsbChwb2ludHMsIHJhbmRvbVNlYXJjaCkge1xuICB2YXIgbiA9IHBvaW50cy5sZW5ndGhcbiAgaWYobiA9PT0gMCkge1xuICAgIHRocm93IG5ldyBFcnJvcihcIk11c3QgaGF2ZSBhdCBsZWFzdCBkKzEgcG9pbnRzXCIpXG4gIH1cbiAgdmFyIGQgPSBwb2ludHNbMF0ubGVuZ3RoXG4gIGlmKG4gPD0gZCkge1xuICAgIHRocm93IG5ldyBFcnJvcihcIk11c3QgaW5wdXQgYXQgbGVhc3QgZCsxIHBvaW50c1wiKVxuICB9XG5cbiAgLy9GSVhNRTogVGhpcyBjb3VsZCBiZSBkZWdlbmVyYXRlLCBidXQgbmVlZCB0byBzZWxlY3QgZCsxIG5vbi1jb3BsYW5hciBwb2ludHMgdG8gYm9vdHN0cmFwIHByb2Nlc3NcbiAgdmFyIGluaXRpYWxTaW1wbGV4ID0gcG9pbnRzLnNsaWNlKDAsIGQrMSlcblxuICAvL01ha2Ugc3VyZSBpbml0aWFsIHNpbXBsZXggaXMgcG9zaXRpdmVseSBvcmllbnRlZFxuICB2YXIgbyA9IG9yaWVudC5hcHBseSh2b2lkIDAsIGluaXRpYWxTaW1wbGV4KVxuICBpZihvID09PSAwKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKFwiSW5wdXQgbm90IGluIGdlbmVyYWwgcG9zaXRpb25cIilcbiAgfVxuICB2YXIgaW5pdGlhbENvb3JkcyA9IG5ldyBBcnJheShkKzEpXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICBpbml0aWFsQ29vcmRzW2ldID0gaVxuICB9XG4gIGlmKG8gPCAwKSB7XG4gICAgaW5pdGlhbENvb3Jkc1swXSA9IDFcbiAgICBpbml0aWFsQ29vcmRzWzFdID0gMFxuICB9XG5cbiAgLy9DcmVhdGUgaW5pdGlhbCB0b3BvbG9naWNhbCBpbmRleCwgZ2x1ZSBwb2ludGVycyB0b2dldGhlciAoa2luZCBvZiBtZXNzeSlcbiAgdmFyIGluaXRpYWxDZWxsID0gbmV3IFNpbXBsZXgoaW5pdGlhbENvb3JkcywgbmV3IEFycmF5KGQrMSksIGZhbHNlKVxuICB2YXIgYm91bmRhcnkgPSBpbml0aWFsQ2VsbC5hZGphY2VudFxuICB2YXIgbGlzdCA9IG5ldyBBcnJheShkKzIpXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgdmVydHMgPSBpbml0aWFsQ29vcmRzLnNsaWNlKClcbiAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICBpZihqID09PSBpKSB7XG4gICAgICAgIHZlcnRzW2pdID0gLTFcbiAgICAgIH1cbiAgICB9XG4gICAgdmFyIHQgPSB2ZXJ0c1swXVxuICAgIHZlcnRzWzBdID0gdmVydHNbMV1cbiAgICB2ZXJ0c1sxXSA9IHRcbiAgICB2YXIgY2VsbCA9IG5ldyBTaW1wbGV4KHZlcnRzLCBuZXcgQXJyYXkoZCsxKSwgdHJ1ZSlcbiAgICBib3VuZGFyeVtpXSA9IGNlbGxcbiAgICBsaXN0W2ldID0gY2VsbFxuICB9XG4gIGxpc3RbZCsxXSA9IGluaXRpYWxDZWxsXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgdmVydHMgPSBib3VuZGFyeVtpXS52ZXJ0aWNlc1xuICAgIHZhciBhZGogPSBib3VuZGFyeVtpXS5hZGphY2VudFxuICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgIHZhciB2ID0gdmVydHNbal1cbiAgICAgIGlmKHYgPCAwKSB7XG4gICAgICAgIGFkaltqXSA9IGluaXRpYWxDZWxsXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBmb3IodmFyIGs9MDsgazw9ZDsgKytrKSB7XG4gICAgICAgIGlmKGJvdW5kYXJ5W2tdLnZlcnRpY2VzLmluZGV4T2YodikgPCAwKSB7XG4gICAgICAgICAgYWRqW2pdID0gYm91bmRhcnlba11cbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgfVxuXG4gIC8vSW5pdGlhbGl6ZSB0cmlhbmdsZXNcbiAgdmFyIHRyaWFuZ2xlcyA9IG5ldyBUcmlhbmd1bGF0aW9uKGQsIGluaXRpYWxTaW1wbGV4LCBsaXN0KVxuXG4gIC8vSW5zZXJ0IHJlbWFpbmluZyBwb2ludHNcbiAgdmFyIHVzZVJhbmRvbSA9ICEhcmFuZG9tU2VhcmNoXG4gIGZvcih2YXIgaT1kKzE7IGk8bjsgKytpKSB7XG4gICAgdHJpYW5nbGVzLmluc2VydChwb2ludHNbaV0sIHVzZVJhbmRvbSlcbiAgfVxuICBcbiAgLy9FeHRyYWN0IGJvdW5kYXJ5IGNlbGxzXG4gIHJldHVybiB0cmlhbmdsZXMuYm91bmRhcnkoKVxufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gZmFzdFR3b1N1bVxuXG5mdW5jdGlvbiBmYXN0VHdvU3VtKGEsIGIsIHJlc3VsdCkge1xuXHR2YXIgeCA9IGEgKyBiXG5cdHZhciBidiA9IHggLSBhXG5cdHZhciBhdiA9IHggLSBidlxuXHR2YXIgYnIgPSBiIC0gYnZcblx0dmFyIGFyID0gYSAtIGF2XG5cdGlmKHJlc3VsdCkge1xuXHRcdHJlc3VsdFswXSA9IGFyICsgYnJcblx0XHRyZXN1bHRbMV0gPSB4XG5cdFx0cmV0dXJuIHJlc3VsdFxuXHR9XG5cdHJldHVybiBbYXIrYnIsIHhdXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxudmFyIHR3b1Byb2R1Y3QgPSByZXF1aXJlKFwidHdvLXByb2R1Y3RcIilcbnZhciB0d29TdW0gPSByZXF1aXJlKFwidHdvLXN1bVwiKVxuXG5tb2R1bGUuZXhwb3J0cyA9IHNjYWxlTGluZWFyRXhwYW5zaW9uXG5cbmZ1bmN0aW9uIHNjYWxlTGluZWFyRXhwYW5zaW9uKGUsIHNjYWxlKSB7XG4gIHZhciBuID0gZS5sZW5ndGhcbiAgaWYobiA9PT0gMSkge1xuICAgIHZhciB0cyA9IHR3b1Byb2R1Y3QoZVswXSwgc2NhbGUpXG4gICAgaWYodHNbMF0pIHtcbiAgICAgIHJldHVybiB0c1xuICAgIH1cbiAgICByZXR1cm4gWyB0c1sxXSBdXG4gIH1cbiAgdmFyIGcgPSBuZXcgQXJyYXkoMiAqIG4pXG4gIHZhciBxID0gWzAuMSwgMC4xXVxuICB2YXIgdCA9IFswLjEsIDAuMV1cbiAgdmFyIGNvdW50ID0gMFxuICB0d29Qcm9kdWN0KGVbMF0sIHNjYWxlLCBxKVxuICBpZihxWzBdKSB7XG4gICAgZ1tjb3VudCsrXSA9IHFbMF1cbiAgfVxuICBmb3IodmFyIGk9MTsgaTxuOyArK2kpIHtcbiAgICB0d29Qcm9kdWN0KGVbaV0sIHNjYWxlLCB0KVxuICAgIHZhciBwcSA9IHFbMV1cbiAgICB0d29TdW0ocHEsIHRbMF0sIHEpXG4gICAgaWYocVswXSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHFbMF1cbiAgICB9XG4gICAgdmFyIGEgPSB0WzFdXG4gICAgdmFyIGIgPSBxWzFdXG4gICAgdmFyIHggPSBhICsgYlxuICAgIHZhciBidiA9IHggLSBhXG4gICAgdmFyIHkgPSBiIC0gYnZcbiAgICBxWzFdID0geFxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICB9XG4gIGlmKHFbMV0pIHtcbiAgICBnW2NvdW50KytdID0gcVsxXVxuICB9XG4gIGlmKGNvdW50ID09PSAwKSB7XG4gICAgZ1tjb3VudCsrXSA9IDAuMFxuICB9XG4gIGcubGVuZ3RoID0gY291bnRcbiAgcmV0dXJuIGdcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IHJvYnVzdFN1YnRyYWN0XG5cbi8vRWFzeSBjYXNlOiBBZGQgdHdvIHNjYWxhcnNcbmZ1bmN0aW9uIHNjYWxhclNjYWxhcihhLCBiKSB7XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIGF2ID0geCAtIGJ2XG4gIHZhciBiciA9IGIgLSBidlxuICB2YXIgYXIgPSBhIC0gYXZcbiAgdmFyIHkgPSBhciArIGJyXG4gIGlmKHkpIHtcbiAgICByZXR1cm4gW3ksIHhdXG4gIH1cbiAgcmV0dXJuIFt4XVxufVxuXG5mdW5jdGlvbiByb2J1c3RTdWJ0cmFjdChlLCBmKSB7XG4gIHZhciBuZSA9IGUubGVuZ3RofDBcbiAgdmFyIG5mID0gZi5sZW5ndGh8MFxuICBpZihuZSA9PT0gMSAmJiBuZiA9PT0gMSkge1xuICAgIHJldHVybiBzY2FsYXJTY2FsYXIoZVswXSwgLWZbMF0pXG4gIH1cbiAgdmFyIG4gPSBuZSArIG5mXG4gIHZhciBnID0gbmV3IEFycmF5KG4pXG4gIHZhciBjb3VudCA9IDBcbiAgdmFyIGVwdHIgPSAwXG4gIHZhciBmcHRyID0gMFxuICB2YXIgYWJzID0gTWF0aC5hYnNcbiAgdmFyIGVpID0gZVtlcHRyXVxuICB2YXIgZWEgPSBhYnMoZWkpXG4gIHZhciBmaSA9IC1mW2ZwdHJdXG4gIHZhciBmYSA9IGFicyhmaSlcbiAgdmFyIGEsIGJcbiAgaWYoZWEgPCBmYSkge1xuICAgIGIgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYiA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgaWYoKGVwdHIgPCBuZSAmJiBlYSA8IGZhKSB8fCAoZnB0ciA+PSBuZikpIHtcbiAgICBhID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGEgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSAtZltmcHRyXVxuICAgICAgZmEgPSBhYnMoZmkpXG4gICAgfVxuICB9XG4gIHZhciB4ID0gYSArIGJcbiAgdmFyIGJ2ID0geCAtIGFcbiAgdmFyIHkgPSBiIC0gYnZcbiAgdmFyIHEwID0geVxuICB2YXIgcTEgPSB4XG4gIHZhciBfeCwgX2J2LCBfYXYsIF9iciwgX2FyXG4gIHdoaWxlKGVwdHIgPCBuZSAmJiBmcHRyIDwgbmYpIHtcbiAgICBpZihlYSA8IGZhKSB7XG4gICAgICBhID0gZWlcbiAgICAgIGVwdHIgKz0gMVxuICAgICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgICBlYSA9IGFicyhlaSlcbiAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgYSA9IGZpXG4gICAgICBmcHRyICs9IDFcbiAgICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgICAgIGZhID0gYWJzKGZpKVxuICAgICAgfVxuICAgIH1cbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gIH1cbiAgd2hpbGUoZXB0ciA8IG5lKSB7XG4gICAgYSA9IGVpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgfVxuICB9XG4gIHdoaWxlKGZwdHIgPCBuZikge1xuICAgIGEgPSBmaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9IFxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgfVxuICB9XG4gIGlmKHEwKSB7XG4gICAgZ1tjb3VudCsrXSA9IHEwXG4gIH1cbiAgaWYocTEpIHtcbiAgICBnW2NvdW50KytdID0gcTFcbiAgfVxuICBpZighY291bnQpIHtcbiAgICBnW2NvdW50KytdID0gMC4wICBcbiAgfVxuICBnLmxlbmd0aCA9IGNvdW50XG4gIHJldHVybiBnXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSBsaW5lYXJFeHBhbnNpb25TdW1cblxuLy9FYXN5IGNhc2U6IEFkZCB0d28gc2NhbGFyc1xuZnVuY3Rpb24gc2NhbGFyU2NhbGFyKGEsIGIpIHtcbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgYXYgPSB4IC0gYnZcbiAgdmFyIGJyID0gYiAtIGJ2XG4gIHZhciBhciA9IGEgLSBhdlxuICB2YXIgeSA9IGFyICsgYnJcbiAgaWYoeSkge1xuICAgIHJldHVybiBbeSwgeF1cbiAgfVxuICByZXR1cm4gW3hdXG59XG5cbmZ1bmN0aW9uIGxpbmVhckV4cGFuc2lvblN1bShlLCBmKSB7XG4gIHZhciBuZSA9IGUubGVuZ3RofDBcbiAgdmFyIG5mID0gZi5sZW5ndGh8MFxuICBpZihuZSA9PT0gMSAmJiBuZiA9PT0gMSkge1xuICAgIHJldHVybiBzY2FsYXJTY2FsYXIoZVswXSwgZlswXSlcbiAgfVxuICB2YXIgbiA9IG5lICsgbmZcbiAgdmFyIGcgPSBuZXcgQXJyYXkobilcbiAgdmFyIGNvdW50ID0gMFxuICB2YXIgZXB0ciA9IDBcbiAgdmFyIGZwdHIgPSAwXG4gIHZhciBhYnMgPSBNYXRoLmFic1xuICB2YXIgZWkgPSBlW2VwdHJdXG4gIHZhciBlYSA9IGFicyhlaSlcbiAgdmFyIGZpID0gZltmcHRyXVxuICB2YXIgZmEgPSBhYnMoZmkpXG4gIHZhciBhLCBiXG4gIGlmKGVhIDwgZmEpIHtcbiAgICBiID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGIgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgaWYoKGVwdHIgPCBuZSAmJiBlYSA8IGZhKSB8fCAoZnB0ciA+PSBuZikpIHtcbiAgICBhID0gZWlcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgICAgZWEgPSBhYnMoZWkpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGEgPSBmaVxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgeSA9IGIgLSBidlxuICB2YXIgcTAgPSB5XG4gIHZhciBxMSA9IHhcbiAgdmFyIF94LCBfYnYsIF9hdiwgX2JyLCBfYXJcbiAgd2hpbGUoZXB0ciA8IG5lICYmIGZwdHIgPCBuZikge1xuICAgIGlmKGVhIDwgZmEpIHtcbiAgICAgIGEgPSBlaVxuICAgICAgZXB0ciArPSAxXG4gICAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICAgIGVhID0gYWJzKGVpKVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICBhID0gZmlcbiAgICAgIGZwdHIgKz0gMVxuICAgICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICAgIGZpID0gZltmcHRyXVxuICAgICAgICBmYSA9IGFicyhmaSlcbiAgICAgIH1cbiAgICB9XG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICB9XG4gIHdoaWxlKGVwdHIgPCBuZSkge1xuICAgIGEgPSBlaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBlcHRyICs9IDFcbiAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgIGVpID0gZVtlcHRyXVxuICAgIH1cbiAgfVxuICB3aGlsZShmcHRyIDwgbmYpIHtcbiAgICBhID0gZmlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfSBcbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGZwdHIgKz0gMVxuICAgIGlmKGZwdHIgPCBuZikge1xuICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgfVxuICB9XG4gIGlmKHEwKSB7XG4gICAgZ1tjb3VudCsrXSA9IHEwXG4gIH1cbiAgaWYocTEpIHtcbiAgICBnW2NvdW50KytdID0gcTFcbiAgfVxuICBpZighY291bnQpIHtcbiAgICBnW2NvdW50KytdID0gMC4wICBcbiAgfVxuICBnLmxlbmd0aCA9IGNvdW50XG4gIHJldHVybiBnXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSB0d29Qcm9kdWN0XG5cbnZhciBTUExJVFRFUiA9ICsoTWF0aC5wb3coMiwgMjcpICsgMS4wKVxuXG5mdW5jdGlvbiB0d29Qcm9kdWN0KGEsIGIsIHJlc3VsdCkge1xuICB2YXIgeCA9IGEgKiBiXG5cbiAgdmFyIGMgPSBTUExJVFRFUiAqIGFcbiAgdmFyIGFiaWcgPSBjIC0gYVxuICB2YXIgYWhpID0gYyAtIGFiaWdcbiAgdmFyIGFsbyA9IGEgLSBhaGlcblxuICB2YXIgZCA9IFNQTElUVEVSICogYlxuICB2YXIgYmJpZyA9IGQgLSBiXG4gIHZhciBiaGkgPSBkIC0gYmJpZ1xuICB2YXIgYmxvID0gYiAtIGJoaVxuXG4gIHZhciBlcnIxID0geCAtIChhaGkgKiBiaGkpXG4gIHZhciBlcnIyID0gZXJyMSAtIChhbG8gKiBiaGkpXG4gIHZhciBlcnIzID0gZXJyMiAtIChhaGkgKiBibG8pXG5cbiAgdmFyIHkgPSBhbG8gKiBibG8gLSBlcnIzXG5cbiAgaWYocmVzdWx0KSB7XG4gICAgcmVzdWx0WzBdID0geVxuICAgIHJlc3VsdFsxXSA9IHhcbiAgICByZXR1cm4gcmVzdWx0XG4gIH1cblxuICByZXR1cm4gWyB5LCB4IF1cbn0iLCJcInVzZSBzdHJpY3RcIlxuXG52YXIgdHdvUHJvZHVjdCA9IHJlcXVpcmUoXCJ0d28tcHJvZHVjdFwiKVxudmFyIHJvYnVzdFN1bSA9IHJlcXVpcmUoXCJyb2J1c3Qtc3VtXCIpXG52YXIgcm9idXN0U2NhbGUgPSByZXF1aXJlKFwicm9idXN0LXNjYWxlXCIpXG52YXIgcm9idXN0U3VidHJhY3QgPSByZXF1aXJlKFwicm9idXN0LXN1YnRyYWN0XCIpXG5cbnZhciBOVU1fRVhQQU5EID0gNVxuXG52YXIgRVBTSUxPTiAgICAgPSAxLjExMDIyMzAyNDYyNTE1NjVlLTE2XG52YXIgRVJSQk9VTkQzICAgPSAoMy4wICsgMTYuMCAqIEVQU0lMT04pICogRVBTSUxPTlxudmFyIEVSUkJPVU5ENCAgID0gKDcuMCArIDU2LjAgKiBFUFNJTE9OKSAqIEVQU0lMT05cblxuZnVuY3Rpb24gY29mYWN0b3IobSwgYykge1xuICB2YXIgcmVzdWx0ID0gbmV3IEFycmF5KG0ubGVuZ3RoLTEpXG4gIGZvcih2YXIgaT0xOyBpPG0ubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgciA9IHJlc3VsdFtpLTFdID0gbmV3IEFycmF5KG0ubGVuZ3RoLTEpXG4gICAgZm9yKHZhciBqPTAsaz0wOyBqPG0ubGVuZ3RoOyArK2opIHtcbiAgICAgIGlmKGogPT09IGMpIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHJbaysrXSA9IG1baV1bal1cbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJlc3VsdFxufVxuXG5mdW5jdGlvbiBtYXRyaXgobikge1xuICB2YXIgcmVzdWx0ID0gbmV3IEFycmF5KG4pXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIHJlc3VsdFtpXSA9IG5ldyBBcnJheShuKVxuICAgIGZvcih2YXIgaj0wOyBqPG47ICsraikge1xuICAgICAgcmVzdWx0W2ldW2pdID0gW1wibVwiLCBqLCBcIltcIiwgKG4taS0xKSwgXCJdXCJdLmpvaW4oXCJcIilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJlc3VsdFxufVxuXG5mdW5jdGlvbiBzaWduKG4pIHtcbiAgaWYobiAmIDEpIHtcbiAgICByZXR1cm4gXCItXCJcbiAgfVxuICByZXR1cm4gXCJcIlxufVxuXG5mdW5jdGlvbiBnZW5lcmF0ZVN1bShleHByKSB7XG4gIGlmKGV4cHIubGVuZ3RoID09PSAxKSB7XG4gICAgcmV0dXJuIGV4cHJbMF1cbiAgfSBlbHNlIGlmKGV4cHIubGVuZ3RoID09PSAyKSB7XG4gICAgcmV0dXJuIFtcInN1bShcIiwgZXhwclswXSwgXCIsXCIsIGV4cHJbMV0sIFwiKVwiXS5qb2luKFwiXCIpXG4gIH0gZWxzZSB7XG4gICAgdmFyIG0gPSBleHByLmxlbmd0aD4+MVxuICAgIHJldHVybiBbXCJzdW0oXCIsIGdlbmVyYXRlU3VtKGV4cHIuc2xpY2UoMCwgbSkpLCBcIixcIiwgZ2VuZXJhdGVTdW0oZXhwci5zbGljZShtKSksIFwiKVwiXS5qb2luKFwiXCIpXG4gIH1cbn1cblxuZnVuY3Rpb24gZGV0ZXJtaW5hbnQobSkge1xuICBpZihtLmxlbmd0aCA9PT0gMikge1xuICAgIHJldHVybiBbW1wic3VtKHByb2QoXCIsIG1bMF1bMF0sIFwiLFwiLCBtWzFdWzFdLCBcIikscHJvZCgtXCIsIG1bMF1bMV0sIFwiLFwiLCBtWzFdWzBdLCBcIikpXCJdLmpvaW4oXCJcIildXG4gIH0gZWxzZSB7XG4gICAgdmFyIGV4cHIgPSBbXVxuICAgIGZvcih2YXIgaT0wOyBpPG0ubGVuZ3RoOyArK2kpIHtcbiAgICAgIGV4cHIucHVzaChbXCJzY2FsZShcIiwgZ2VuZXJhdGVTdW0oZGV0ZXJtaW5hbnQoY29mYWN0b3IobSwgaSkpKSwgXCIsXCIsIHNpZ24oaSksIG1bMF1baV0sIFwiKVwiXS5qb2luKFwiXCIpKVxuICAgIH1cbiAgICByZXR1cm4gZXhwclxuICB9XG59XG5cbmZ1bmN0aW9uIG9yaWVudGF0aW9uKG4pIHtcbiAgdmFyIHBvcyA9IFtdXG4gIHZhciBuZWcgPSBbXVxuICB2YXIgbSA9IG1hdHJpeChuKVxuICB2YXIgYXJncyA9IFtdXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIGlmKChpJjEpPT09MCkge1xuICAgICAgcG9zLnB1c2guYXBwbHkocG9zLCBkZXRlcm1pbmFudChjb2ZhY3RvcihtLCBpKSkpXG4gICAgfSBlbHNlIHtcbiAgICAgIG5lZy5wdXNoLmFwcGx5KG5lZywgZGV0ZXJtaW5hbnQoY29mYWN0b3IobSwgaSkpKVxuICAgIH1cbiAgICBhcmdzLnB1c2goXCJtXCIgKyBpKVxuICB9XG4gIHZhciBwb3NFeHByID0gZ2VuZXJhdGVTdW0ocG9zKVxuICB2YXIgbmVnRXhwciA9IGdlbmVyYXRlU3VtKG5lZylcbiAgdmFyIGZ1bmNOYW1lID0gXCJvcmllbnRhdGlvblwiICsgbiArIFwiRXhhY3RcIlxuICB2YXIgY29kZSA9IFtcImZ1bmN0aW9uIFwiLCBmdW5jTmFtZSwgXCIoXCIsIGFyZ3Muam9pbigpLCBcIil7dmFyIHA9XCIsIHBvc0V4cHIsIFwiLG49XCIsIG5lZ0V4cHIsIFwiLGQ9c3ViKHAsbik7XFxcbnJldHVybiBkW2QubGVuZ3RoLTFdO307cmV0dXJuIFwiLCBmdW5jTmFtZV0uam9pbihcIlwiKVxuICB2YXIgcHJvYyA9IG5ldyBGdW5jdGlvbihcInN1bVwiLCBcInByb2RcIiwgXCJzY2FsZVwiLCBcInN1YlwiLCBjb2RlKVxuICByZXR1cm4gcHJvYyhyb2J1c3RTdW0sIHR3b1Byb2R1Y3QsIHJvYnVzdFNjYWxlLCByb2J1c3RTdWJ0cmFjdClcbn1cblxudmFyIG9yaWVudGF0aW9uM0V4YWN0ID0gb3JpZW50YXRpb24oMylcbnZhciBvcmllbnRhdGlvbjRFeGFjdCA9IG9yaWVudGF0aW9uKDQpXG5cbnZhciBDQUNIRUQgPSBbXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMCgpIHsgcmV0dXJuIDAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24xKCkgeyByZXR1cm4gMCB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjIoYSwgYikgeyBcbiAgICByZXR1cm4gYlswXSAtIGFbMF1cbiAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24zKGEsIGIsIGMpIHtcbiAgICB2YXIgbCA9IChhWzFdIC0gY1sxXSkgKiAoYlswXSAtIGNbMF0pXG4gICAgdmFyIHIgPSAoYVswXSAtIGNbMF0pICogKGJbMV0gLSBjWzFdKVxuICAgIHZhciBkZXQgPSBsIC0gclxuICAgIHZhciBzXG4gICAgaWYobCA+IDApIHtcbiAgICAgIGlmKHIgPD0gMCkge1xuICAgICAgICByZXR1cm4gZGV0XG4gICAgICB9IGVsc2Uge1xuICAgICAgICBzID0gbCArIHJcbiAgICAgIH1cbiAgICB9IGVsc2UgaWYobCA8IDApIHtcbiAgICAgIGlmKHIgPj0gMCkge1xuICAgICAgICByZXR1cm4gZGV0XG4gICAgICB9IGVsc2Uge1xuICAgICAgICBzID0gLShsICsgcilcbiAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgcmV0dXJuIGRldFxuICAgIH1cbiAgICB2YXIgdG9sID0gRVJSQk9VTkQzICogc1xuICAgIGlmKGRldCA+PSB0b2wgfHwgZGV0IDw9IC10b2wpIHtcbiAgICAgIHJldHVybiBkZXRcbiAgICB9XG4gICAgcmV0dXJuIG9yaWVudGF0aW9uM0V4YWN0KGEsIGIsIGMpXG4gIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uNChhLGIsYyxkKSB7XG4gICAgdmFyIGFkeCA9IGFbMF0gLSBkWzBdXG4gICAgdmFyIGJkeCA9IGJbMF0gLSBkWzBdXG4gICAgdmFyIGNkeCA9IGNbMF0gLSBkWzBdXG4gICAgdmFyIGFkeSA9IGFbMV0gLSBkWzFdXG4gICAgdmFyIGJkeSA9IGJbMV0gLSBkWzFdXG4gICAgdmFyIGNkeSA9IGNbMV0gLSBkWzFdXG4gICAgdmFyIGFkeiA9IGFbMl0gLSBkWzJdXG4gICAgdmFyIGJkeiA9IGJbMl0gLSBkWzJdXG4gICAgdmFyIGNkeiA9IGNbMl0gLSBkWzJdXG4gICAgdmFyIGJkeGNkeSA9IGJkeCAqIGNkeVxuICAgIHZhciBjZHhiZHkgPSBjZHggKiBiZHlcbiAgICB2YXIgY2R4YWR5ID0gY2R4ICogYWR5XG4gICAgdmFyIGFkeGNkeSA9IGFkeCAqIGNkeVxuICAgIHZhciBhZHhiZHkgPSBhZHggKiBiZHlcbiAgICB2YXIgYmR4YWR5ID0gYmR4ICogYWR5XG4gICAgdmFyIGRldCA9IGFkeiAqIChiZHhjZHkgLSBjZHhiZHkpIFxuICAgICAgICAgICAgKyBiZHogKiAoY2R4YWR5IC0gYWR4Y2R5KVxuICAgICAgICAgICAgKyBjZHogKiAoYWR4YmR5IC0gYmR4YWR5KVxuICAgIHZhciBwZXJtYW5lbnQgPSAoTWF0aC5hYnMoYmR4Y2R5KSArIE1hdGguYWJzKGNkeGJkeSkpICogTWF0aC5hYnMoYWR6KVxuICAgICAgICAgICAgICAgICAgKyAoTWF0aC5hYnMoY2R4YWR5KSArIE1hdGguYWJzKGFkeGNkeSkpICogTWF0aC5hYnMoYmR6KVxuICAgICAgICAgICAgICAgICAgKyAoTWF0aC5hYnMoYWR4YmR5KSArIE1hdGguYWJzKGJkeGFkeSkpICogTWF0aC5hYnMoY2R6KVxuICAgIHZhciB0b2wgPSBFUlJCT1VORDQgKiBwZXJtYW5lbnRcbiAgICBpZiAoKGRldCA+IHRvbCkgfHwgKC1kZXQgPiB0b2wpKSB7XG4gICAgICByZXR1cm4gZGV0XG4gICAgfVxuICAgIHJldHVybiBvcmllbnRhdGlvbjRFeGFjdChhLGIsYyxkKVxuICB9XG5dXG5cbmZ1bmN0aW9uIHNsb3dPcmllbnQoYXJncykge1xuICB2YXIgcHJvYyA9IENBQ0hFRFthcmdzLmxlbmd0aF1cbiAgaWYoIXByb2MpIHtcbiAgICBwcm9jID0gQ0FDSEVEW2FyZ3MubGVuZ3RoXSA9IG9yaWVudGF0aW9uKGFyZ3MubGVuZ3RoKVxuICB9XG4gIHJldHVybiBwcm9jLmFwcGx5KHVuZGVmaW5lZCwgYXJncylcbn1cblxuZnVuY3Rpb24gZ2VuZXJhdGVPcmllbnRhdGlvblByb2MoKSB7XG4gIHdoaWxlKENBQ0hFRC5sZW5ndGggPD0gTlVNX0VYUEFORCkge1xuICAgIENBQ0hFRC5wdXNoKG9yaWVudGF0aW9uKENBQ0hFRC5sZW5ndGgpKVxuICB9XG4gIHZhciBhcmdzID0gW11cbiAgdmFyIHByb2NBcmdzID0gW1wic2xvd1wiXVxuICBmb3IodmFyIGk9MDsgaTw9TlVNX0VYUEFORDsgKytpKSB7XG4gICAgYXJncy5wdXNoKFwiYVwiICsgaSlcbiAgICBwcm9jQXJncy5wdXNoKFwib1wiICsgaSlcbiAgfVxuICB2YXIgY29kZSA9IFtcbiAgICBcImZ1bmN0aW9uIGdldE9yaWVudGF0aW9uKFwiLCBhcmdzLmpvaW4oKSwgXCIpe3N3aXRjaChhcmd1bWVudHMubGVuZ3RoKXtjYXNlIDA6Y2FzZSAxOnJldHVybiAwO1wiXG4gIF1cbiAgZm9yKHZhciBpPTI7IGk8PU5VTV9FWFBBTkQ7ICsraSkge1xuICAgIGNvZGUucHVzaChcImNhc2UgXCIsIGksIFwiOnJldHVybiBvXCIsIGksIFwiKFwiLCBhcmdzLnNsaWNlKDAsIGkpLmpvaW4oKSwgXCIpO1wiKVxuICB9XG4gIGNvZGUucHVzaChcIn12YXIgcz1uZXcgQXJyYXkoYXJndW1lbnRzLmxlbmd0aCk7Zm9yKHZhciBpPTA7aTxhcmd1bWVudHMubGVuZ3RoOysraSl7c1tpXT1hcmd1bWVudHNbaV19O3JldHVybiBzbG93KHMpO31yZXR1cm4gZ2V0T3JpZW50YXRpb25cIilcbiAgcHJvY0FyZ3MucHVzaChjb2RlLmpvaW4oXCJcIikpXG5cbiAgdmFyIHByb2MgPSBGdW5jdGlvbi5hcHBseSh1bmRlZmluZWQsIHByb2NBcmdzKVxuICBtb2R1bGUuZXhwb3J0cyA9IHByb2MuYXBwbHkodW5kZWZpbmVkLCBbc2xvd09yaWVudF0uY29uY2F0KENBQ0hFRCkpXG4gIGZvcih2YXIgaT0wOyBpPD1OVU1fRVhQQU5EOyArK2kpIHtcbiAgICBtb2R1bGUuZXhwb3J0c1tpXSA9IENBQ0hFRFtpXVxuICB9XG59XG5cbmdlbmVyYXRlT3JpZW50YXRpb25Qcm9jKCkiLCIvKipcbiAqIEJpdCB0d2lkZGxpbmcgaGFja3MgZm9yIEphdmFTY3JpcHQuXG4gKlxuICogQXV0aG9yOiBNaWtvbGEgTHlzZW5rb1xuICpcbiAqIFBvcnRlZCBmcm9tIFN0YW5mb3JkIGJpdCB0d2lkZGxpbmcgaGFjayBsaWJyYXJ5OlxuICogICAgaHR0cDovL2dyYXBoaWNzLnN0YW5mb3JkLmVkdS9+c2VhbmRlci9iaXRoYWNrcy5odG1sXG4gKi9cblxuXCJ1c2Ugc3RyaWN0XCI7IFwidXNlIHJlc3RyaWN0XCI7XG5cbi8vTnVtYmVyIG9mIGJpdHMgaW4gYW4gaW50ZWdlclxudmFyIElOVF9CSVRTID0gMzI7XG5cbi8vQ29uc3RhbnRzXG5leHBvcnRzLklOVF9CSVRTICA9IElOVF9CSVRTO1xuZXhwb3J0cy5JTlRfTUFYICAgPSAgMHg3ZmZmZmZmZjtcbmV4cG9ydHMuSU5UX01JTiAgID0gLTE8PChJTlRfQklUUy0xKTtcblxuLy9SZXR1cm5zIC0xLCAwLCArMSBkZXBlbmRpbmcgb24gc2lnbiBvZiB4XG5leHBvcnRzLnNpZ24gPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAodiA+IDApIC0gKHYgPCAwKTtcbn1cblxuLy9Db21wdXRlcyBhYnNvbHV0ZSB2YWx1ZSBvZiBpbnRlZ2VyXG5leHBvcnRzLmFicyA9IGZ1bmN0aW9uKHYpIHtcbiAgdmFyIG1hc2sgPSB2ID4+IChJTlRfQklUUy0xKTtcbiAgcmV0dXJuICh2IF4gbWFzaykgLSBtYXNrO1xufVxuXG4vL0NvbXB1dGVzIG1pbmltdW0gb2YgaW50ZWdlcnMgeCBhbmQgeVxuZXhwb3J0cy5taW4gPSBmdW5jdGlvbih4LCB5KSB7XG4gIHJldHVybiB5IF4gKCh4IF4geSkgJiAtKHggPCB5KSk7XG59XG5cbi8vQ29tcHV0ZXMgbWF4aW11bSBvZiBpbnRlZ2VycyB4IGFuZCB5XG5leHBvcnRzLm1heCA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgcmV0dXJuIHggXiAoKHggXiB5KSAmIC0oeCA8IHkpKTtcbn1cblxuLy9DaGVja3MgaWYgYSBudW1iZXIgaXMgYSBwb3dlciBvZiB0d29cbmV4cG9ydHMuaXNQb3cyID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gISh2ICYgKHYtMSkpICYmICghIXYpO1xufVxuXG4vL0NvbXB1dGVzIGxvZyBiYXNlIDIgb2YgdlxuZXhwb3J0cy5sb2cyID0gZnVuY3Rpb24odikge1xuICB2YXIgciwgc2hpZnQ7XG4gIHIgPSAgICAgKHYgPiAweEZGRkYpIDw8IDQ7IHYgPj4+PSByO1xuICBzaGlmdCA9ICh2ID4gMHhGRiAgKSA8PCAzOyB2ID4+Pj0gc2hpZnQ7IHIgfD0gc2hpZnQ7XG4gIHNoaWZ0ID0gKHYgPiAweEYgICApIDw8IDI7IHYgPj4+PSBzaGlmdDsgciB8PSBzaGlmdDtcbiAgc2hpZnQgPSAodiA+IDB4MyAgICkgPDwgMTsgdiA+Pj49IHNoaWZ0OyByIHw9IHNoaWZ0O1xuICByZXR1cm4gciB8ICh2ID4+IDEpO1xufVxuXG4vL0NvbXB1dGVzIGxvZyBiYXNlIDEwIG9mIHZcbmV4cG9ydHMubG9nMTAgPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAgKHYgPj0gMTAwMDAwMDAwMCkgPyA5IDogKHYgPj0gMTAwMDAwMDAwKSA/IDggOiAodiA+PSAxMDAwMDAwMCkgPyA3IDpcbiAgICAgICAgICAodiA+PSAxMDAwMDAwKSA/IDYgOiAodiA+PSAxMDAwMDApID8gNSA6ICh2ID49IDEwMDAwKSA/IDQgOlxuICAgICAgICAgICh2ID49IDEwMDApID8gMyA6ICh2ID49IDEwMCkgPyAyIDogKHYgPj0gMTApID8gMSA6IDA7XG59XG5cbi8vQ291bnRzIG51bWJlciBvZiBiaXRzXG5leHBvcnRzLnBvcENvdW50ID0gZnVuY3Rpb24odikge1xuICB2ID0gdiAtICgodiA+Pj4gMSkgJiAweDU1NTU1NTU1KTtcbiAgdiA9ICh2ICYgMHgzMzMzMzMzMykgKyAoKHYgPj4+IDIpICYgMHgzMzMzMzMzMyk7XG4gIHJldHVybiAoKHYgKyAodiA+Pj4gNCkgJiAweEYwRjBGMEYpICogMHgxMDEwMTAxKSA+Pj4gMjQ7XG59XG5cbi8vQ291bnRzIG51bWJlciBvZiB0cmFpbGluZyB6ZXJvc1xuZnVuY3Rpb24gY291bnRUcmFpbGluZ1plcm9zKHYpIHtcbiAgdmFyIGMgPSAzMjtcbiAgdiAmPSAtdjtcbiAgaWYgKHYpIGMtLTtcbiAgaWYgKHYgJiAweDAwMDBGRkZGKSBjIC09IDE2O1xuICBpZiAodiAmIDB4MDBGRjAwRkYpIGMgLT0gODtcbiAgaWYgKHYgJiAweDBGMEYwRjBGKSBjIC09IDQ7XG4gIGlmICh2ICYgMHgzMzMzMzMzMykgYyAtPSAyO1xuICBpZiAodiAmIDB4NTU1NTU1NTUpIGMgLT0gMTtcbiAgcmV0dXJuIGM7XG59XG5leHBvcnRzLmNvdW50VHJhaWxpbmdaZXJvcyA9IGNvdW50VHJhaWxpbmdaZXJvcztcblxuLy9Sb3VuZHMgdG8gbmV4dCBwb3dlciBvZiAyXG5leHBvcnRzLm5leHRQb3cyID0gZnVuY3Rpb24odikge1xuICB2ICs9IHYgPT09IDA7XG4gIC0tdjtcbiAgdiB8PSB2ID4+PiAxO1xuICB2IHw9IHYgPj4+IDI7XG4gIHYgfD0gdiA+Pj4gNDtcbiAgdiB8PSB2ID4+PiA4O1xuICB2IHw9IHYgPj4+IDE2O1xuICByZXR1cm4gdiArIDE7XG59XG5cbi8vUm91bmRzIGRvd24gdG8gcHJldmlvdXMgcG93ZXIgb2YgMlxuZXhwb3J0cy5wcmV2UG93MiA9IGZ1bmN0aW9uKHYpIHtcbiAgdiB8PSB2ID4+PiAxO1xuICB2IHw9IHYgPj4+IDI7XG4gIHYgfD0gdiA+Pj4gNDtcbiAgdiB8PSB2ID4+PiA4O1xuICB2IHw9IHYgPj4+IDE2O1xuICByZXR1cm4gdiAtICh2Pj4+MSk7XG59XG5cbi8vQ29tcHV0ZXMgcGFyaXR5IG9mIHdvcmRcbmV4cG9ydHMucGFyaXR5ID0gZnVuY3Rpb24odikge1xuICB2IF49IHYgPj4+IDE2O1xuICB2IF49IHYgPj4+IDg7XG4gIHYgXj0gdiA+Pj4gNDtcbiAgdiAmPSAweGY7XG4gIHJldHVybiAoMHg2OTk2ID4+PiB2KSAmIDE7XG59XG5cbnZhciBSRVZFUlNFX1RBQkxFID0gbmV3IEFycmF5KDI1Nik7XG5cbihmdW5jdGlvbih0YWIpIHtcbiAgZm9yKHZhciBpPTA7IGk8MjU2OyArK2kpIHtcbiAgICB2YXIgdiA9IGksIHIgPSBpLCBzID0gNztcbiAgICBmb3IgKHYgPj4+PSAxOyB2OyB2ID4+Pj0gMSkge1xuICAgICAgciA8PD0gMTtcbiAgICAgIHIgfD0gdiAmIDE7XG4gICAgICAtLXM7XG4gICAgfVxuICAgIHRhYltpXSA9IChyIDw8IHMpICYgMHhmZjtcbiAgfVxufSkoUkVWRVJTRV9UQUJMRSk7XG5cbi8vUmV2ZXJzZSBiaXRzIGluIGEgMzIgYml0IHdvcmRcbmV4cG9ydHMucmV2ZXJzZSA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICAoUkVWRVJTRV9UQUJMRVsgdiAgICAgICAgICYgMHhmZl0gPDwgMjQpIHxcbiAgICAgICAgICAoUkVWRVJTRV9UQUJMRVsodiA+Pj4gOCkgICYgMHhmZl0gPDwgMTYpIHxcbiAgICAgICAgICAoUkVWRVJTRV9UQUJMRVsodiA+Pj4gMTYpICYgMHhmZl0gPDwgOCkgIHxcbiAgICAgICAgICAgUkVWRVJTRV9UQUJMRVsodiA+Pj4gMjQpICYgMHhmZl07XG59XG5cbi8vSW50ZXJsZWF2ZSBiaXRzIG9mIDIgY29vcmRpbmF0ZXMgd2l0aCAxNiBiaXRzLiAgVXNlZnVsIGZvciBmYXN0IHF1YWR0cmVlIGNvZGVzXG5leHBvcnRzLmludGVybGVhdmUyID0gZnVuY3Rpb24oeCwgeSkge1xuICB4ICY9IDB4RkZGRjtcbiAgeCA9ICh4IHwgKHggPDwgOCkpICYgMHgwMEZGMDBGRjtcbiAgeCA9ICh4IHwgKHggPDwgNCkpICYgMHgwRjBGMEYwRjtcbiAgeCA9ICh4IHwgKHggPDwgMikpICYgMHgzMzMzMzMzMztcbiAgeCA9ICh4IHwgKHggPDwgMSkpICYgMHg1NTU1NTU1NTtcblxuICB5ICY9IDB4RkZGRjtcbiAgeSA9ICh5IHwgKHkgPDwgOCkpICYgMHgwMEZGMDBGRjtcbiAgeSA9ICh5IHwgKHkgPDwgNCkpICYgMHgwRjBGMEYwRjtcbiAgeSA9ICh5IHwgKHkgPDwgMikpICYgMHgzMzMzMzMzMztcbiAgeSA9ICh5IHwgKHkgPDwgMSkpICYgMHg1NTU1NTU1NTtcblxuICByZXR1cm4geCB8ICh5IDw8IDEpO1xufVxuXG4vL0V4dHJhY3RzIHRoZSBudGggaW50ZXJsZWF2ZWQgY29tcG9uZW50XG5leHBvcnRzLmRlaW50ZXJsZWF2ZTIgPSBmdW5jdGlvbih2LCBuKSB7XG4gIHYgPSAodiA+Pj4gbikgJiAweDU1NTU1NTU1O1xuICB2ID0gKHYgfCAodiA+Pj4gMSkpICAmIDB4MzMzMzMzMzM7XG4gIHYgPSAodiB8ICh2ID4+PiAyKSkgICYgMHgwRjBGMEYwRjtcbiAgdiA9ICh2IHwgKHYgPj4+IDQpKSAgJiAweDAwRkYwMEZGO1xuICB2ID0gKHYgfCAodiA+Pj4gMTYpKSAmIDB4MDAwRkZGRjtcbiAgcmV0dXJuICh2IDw8IDE2KSA+PiAxNjtcbn1cblxuXG4vL0ludGVybGVhdmUgYml0cyBvZiAzIGNvb3JkaW5hdGVzLCBlYWNoIHdpdGggMTAgYml0cy4gIFVzZWZ1bCBmb3IgZmFzdCBvY3RyZWUgY29kZXNcbmV4cG9ydHMuaW50ZXJsZWF2ZTMgPSBmdW5jdGlvbih4LCB5LCB6KSB7XG4gIHggJj0gMHgzRkY7XG4gIHggID0gKHggfCAoeDw8MTYpKSAmIDQyNzgxOTAzMzU7XG4gIHggID0gKHggfCAoeDw8OCkpICAmIDI1MTcxOTY5NTtcbiAgeCAgPSAoeCB8ICh4PDw0KSkgICYgMzI3MjM1NjAzNTtcbiAgeCAgPSAoeCB8ICh4PDwyKSkgICYgMTIyNzEzMzUxMztcblxuICB5ICY9IDB4M0ZGO1xuICB5ICA9ICh5IHwgKHk8PDE2KSkgJiA0Mjc4MTkwMzM1O1xuICB5ICA9ICh5IHwgKHk8PDgpKSAgJiAyNTE3MTk2OTU7XG4gIHkgID0gKHkgfCAoeTw8NCkpICAmIDMyNzIzNTYwMzU7XG4gIHkgID0gKHkgfCAoeTw8MikpICAmIDEyMjcxMzM1MTM7XG4gIHggfD0gKHkgPDwgMSk7XG4gIFxuICB6ICY9IDB4M0ZGO1xuICB6ICA9ICh6IHwgKHo8PDE2KSkgJiA0Mjc4MTkwMzM1O1xuICB6ICA9ICh6IHwgKHo8PDgpKSAgJiAyNTE3MTk2OTU7XG4gIHogID0gKHogfCAoejw8NCkpICAmIDMyNzIzNTYwMzU7XG4gIHogID0gKHogfCAoejw8MikpICAmIDEyMjcxMzM1MTM7XG4gIFxuICByZXR1cm4geCB8ICh6IDw8IDIpO1xufVxuXG4vL0V4dHJhY3RzIG50aCBpbnRlcmxlYXZlZCBjb21wb25lbnQgb2YgYSAzLXR1cGxlXG5leHBvcnRzLmRlaW50ZXJsZWF2ZTMgPSBmdW5jdGlvbih2LCBuKSB7XG4gIHYgPSAodiA+Pj4gbikgICAgICAgJiAxMjI3MTMzNTEzO1xuICB2ID0gKHYgfCAodj4+PjIpKSAgICYgMzI3MjM1NjAzNTtcbiAgdiA9ICh2IHwgKHY+Pj40KSkgICAmIDI1MTcxOTY5NTtcbiAgdiA9ICh2IHwgKHY+Pj44KSkgICAmIDQyNzgxOTAzMzU7XG4gIHYgPSAodiB8ICh2Pj4+MTYpKSAgJiAweDNGRjtcbiAgcmV0dXJuICh2PDwyMik+PjIyO1xufVxuXG4vL0NvbXB1dGVzIG5leHQgY29tYmluYXRpb24gaW4gY29sZXhpY29ncmFwaGljIG9yZGVyICh0aGlzIGlzIG1pc3Rha2VubHkgY2FsbGVkIG5leHRQZXJtdXRhdGlvbiBvbiB0aGUgYml0IHR3aWRkbGluZyBoYWNrcyBwYWdlKVxuZXhwb3J0cy5uZXh0Q29tYmluYXRpb24gPSBmdW5jdGlvbih2KSB7XG4gIHZhciB0ID0gdiB8ICh2IC0gMSk7XG4gIHJldHVybiAodCArIDEpIHwgKCgofnQgJiAtfnQpIC0gMSkgPj4+IChjb3VudFRyYWlsaW5nWmVyb3ModikgKyAxKSk7XG59XG5cbiIsIlwidXNlIHN0cmljdFwiOyBcInVzZSByZXN0cmljdFwiO1xuXG5tb2R1bGUuZXhwb3J0cyA9IFVuaW9uRmluZDtcblxuZnVuY3Rpb24gVW5pb25GaW5kKGNvdW50KSB7XG4gIHRoaXMucm9vdHMgPSBuZXcgQXJyYXkoY291bnQpO1xuICB0aGlzLnJhbmtzID0gbmV3IEFycmF5KGNvdW50KTtcbiAgXG4gIGZvcih2YXIgaT0wOyBpPGNvdW50OyArK2kpIHtcbiAgICB0aGlzLnJvb3RzW2ldID0gaTtcbiAgICB0aGlzLnJhbmtzW2ldID0gMDtcbiAgfVxufVxuXG52YXIgcHJvdG8gPSBVbmlvbkZpbmQucHJvdG90eXBlXG5cbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShwcm90bywgXCJsZW5ndGhcIiwge1xuICBcImdldFwiOiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5yb290cy5sZW5ndGhcbiAgfVxufSlcblxucHJvdG8ubWFrZVNldCA9IGZ1bmN0aW9uKCkge1xuICB2YXIgbiA9IHRoaXMucm9vdHMubGVuZ3RoO1xuICB0aGlzLnJvb3RzLnB1c2gobik7XG4gIHRoaXMucmFua3MucHVzaCgwKTtcbiAgcmV0dXJuIG47XG59XG5cbnByb3RvLmZpbmQgPSBmdW5jdGlvbih4KSB7XG4gIHZhciByb290cyA9IHRoaXMucm9vdHM7XG4gIHdoaWxlKHJvb3RzW3hdICE9PSB4KSB7XG4gICAgdmFyIHkgPSByb290c1t4XTtcbiAgICByb290c1t4XSA9IHJvb3RzW3ldO1xuICAgIHggPSB5O1xuICB9XG4gIHJldHVybiB4O1xufVxuXG5wcm90by5saW5rID0gZnVuY3Rpb24oeCwgeSkge1xuICB2YXIgeHIgPSB0aGlzLmZpbmQoeClcbiAgICAsIHlyID0gdGhpcy5maW5kKHkpO1xuICBpZih4ciA9PT0geXIpIHtcbiAgICByZXR1cm47XG4gIH1cbiAgdmFyIHJhbmtzID0gdGhpcy5yYW5rc1xuICAgICwgcm9vdHMgPSB0aGlzLnJvb3RzXG4gICAgLCB4ZCAgICA9IHJhbmtzW3hyXVxuICAgICwgeWQgICAgPSByYW5rc1t5cl07XG4gIGlmKHhkIDwgeWQpIHtcbiAgICByb290c1t4cl0gPSB5cjtcbiAgfSBlbHNlIGlmKHlkIDwgeGQpIHtcbiAgICByb290c1t5cl0gPSB4cjtcbiAgfSBlbHNlIHtcbiAgICByb290c1t5cl0gPSB4cjtcbiAgICArK3JhbmtzW3hyXTtcbiAgfVxufSIsIlwidXNlIHN0cmljdFwiOyBcInVzZSByZXN0cmljdFwiO1xuXG52YXIgYml0cyAgICAgID0gcmVxdWlyZShcImJpdC10d2lkZGxlXCIpXG4gICwgVW5pb25GaW5kID0gcmVxdWlyZShcInVuaW9uLWZpbmRcIilcblxuLy9SZXR1cm5zIHRoZSBkaW1lbnNpb24gb2YgYSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGRpbWVuc2lvbihjZWxscykge1xuICB2YXIgZCA9IDBcbiAgICAsIG1heCA9IE1hdGgubWF4XG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIGQgPSBtYXgoZCwgY2VsbHNbaV0ubGVuZ3RoKVxuICB9XG4gIHJldHVybiBkLTFcbn1cbmV4cG9ydHMuZGltZW5zaW9uID0gZGltZW5zaW9uXG5cbi8vQ291bnRzIHRoZSBudW1iZXIgb2YgdmVydGljZXMgaW4gZmFjZXNcbmZ1bmN0aW9uIGNvdW50VmVydGljZXMoY2VsbHMpIHtcbiAgdmFyIHZjID0gLTFcbiAgICAsIG1heCA9IE1hdGgubWF4XG4gIGZvcih2YXIgaT0wLCBpbD1jZWxscy5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MCwgamw9Yy5sZW5ndGg7IGo8amw7ICsraikge1xuICAgICAgdmMgPSBtYXgodmMsIGNbal0pXG4gICAgfVxuICB9XG4gIHJldHVybiB2YysxXG59XG5leHBvcnRzLmNvdW50VmVydGljZXMgPSBjb3VudFZlcnRpY2VzXG5cbi8vUmV0dXJucyBhIGRlZXAgY29weSBvZiBjZWxsc1xuZnVuY3Rpb24gY2xvbmVDZWxscyhjZWxscykge1xuICB2YXIgbmNlbGxzID0gbmV3IEFycmF5KGNlbGxzLmxlbmd0aClcbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgbmNlbGxzW2ldID0gY2VsbHNbaV0uc2xpY2UoMClcbiAgfVxuICByZXR1cm4gbmNlbGxzXG59XG5leHBvcnRzLmNsb25lQ2VsbHMgPSBjbG9uZUNlbGxzXG5cbi8vUmFua3MgYSBwYWlyIG9mIGNlbGxzIHVwIHRvIHBlcm11dGF0aW9uXG5mdW5jdGlvbiBjb21wYXJlQ2VsbHMoYSwgYikge1xuICB2YXIgbiA9IGEubGVuZ3RoXG4gICAgLCB0ID0gYS5sZW5ndGggLSBiLmxlbmd0aFxuICAgICwgbWluID0gTWF0aC5taW5cbiAgaWYodCkge1xuICAgIHJldHVybiB0XG4gIH1cbiAgc3dpdGNoKG4pIHtcbiAgICBjYXNlIDA6XG4gICAgICByZXR1cm4gMDtcbiAgICBjYXNlIDE6XG4gICAgICByZXR1cm4gYVswXSAtIGJbMF07XG4gICAgY2FzZSAyOlxuICAgICAgdmFyIGQgPSBhWzBdK2FbMV0tYlswXS1iWzFdXG4gICAgICBpZihkKSB7XG4gICAgICAgIHJldHVybiBkXG4gICAgICB9XG4gICAgICByZXR1cm4gbWluKGFbMF0sYVsxXSkgLSBtaW4oYlswXSxiWzFdKVxuICAgIGNhc2UgMzpcbiAgICAgIHZhciBsMSA9IGFbMF0rYVsxXVxuICAgICAgICAsIG0xID0gYlswXStiWzFdXG4gICAgICBkID0gbDErYVsyXSAtIChtMStiWzJdKVxuICAgICAgaWYoZCkge1xuICAgICAgICByZXR1cm4gZFxuICAgICAgfVxuICAgICAgdmFyIGwwID0gbWluKGFbMF0sIGFbMV0pXG4gICAgICAgICwgbTAgPSBtaW4oYlswXSwgYlsxXSlcbiAgICAgICAgLCBkICA9IG1pbihsMCwgYVsyXSkgLSBtaW4obTAsIGJbMl0pXG4gICAgICBpZihkKSB7XG4gICAgICAgIHJldHVybiBkXG4gICAgICB9XG4gICAgICByZXR1cm4gbWluKGwwK2FbMl0sIGwxKSAtIG1pbihtMCtiWzJdLCBtMSlcbiAgICBcbiAgICAvL1RPRE86IE1heWJlIG9wdGltaXplIG49NCBhcyB3ZWxsP1xuICAgIFxuICAgIGRlZmF1bHQ6XG4gICAgICB2YXIgYXMgPSBhLnNsaWNlKDApXG4gICAgICBhcy5zb3J0KClcbiAgICAgIHZhciBicyA9IGIuc2xpY2UoMClcbiAgICAgIGJzLnNvcnQoKVxuICAgICAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgICAgIHQgPSBhc1tpXSAtIGJzW2ldXG4gICAgICAgIGlmKHQpIHtcbiAgICAgICAgICByZXR1cm4gdFxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXR1cm4gMFxuICB9XG59XG5leHBvcnRzLmNvbXBhcmVDZWxscyA9IGNvbXBhcmVDZWxsc1xuXG5mdW5jdGlvbiBjb21wYXJlWmlwcGVkKGEsIGIpIHtcbiAgcmV0dXJuIGNvbXBhcmVDZWxscyhhWzBdLCBiWzBdKVxufVxuXG4vL1B1dHMgYSBjZWxsIGNvbXBsZXggaW50byBub3JtYWwgb3JkZXIgZm9yIHRoZSBwdXJwb3NlcyBvZiBmaW5kQ2VsbCBxdWVyaWVzXG5mdW5jdGlvbiBub3JtYWxpemUoY2VsbHMsIGF0dHIpIHtcbiAgaWYoYXR0cikge1xuICAgIHZhciBsZW4gPSBjZWxscy5sZW5ndGhcbiAgICB2YXIgemlwcGVkID0gbmV3IEFycmF5KGxlbilcbiAgICBmb3IodmFyIGk9MDsgaTxsZW47ICsraSkge1xuICAgICAgemlwcGVkW2ldID0gW2NlbGxzW2ldLCBhdHRyW2ldXVxuICAgIH1cbiAgICB6aXBwZWQuc29ydChjb21wYXJlWmlwcGVkKVxuICAgIGZvcih2YXIgaT0wOyBpPGxlbjsgKytpKSB7XG4gICAgICBjZWxsc1tpXSA9IHppcHBlZFtpXVswXVxuICAgICAgYXR0cltpXSA9IHppcHBlZFtpXVsxXVxuICAgIH1cbiAgICByZXR1cm4gY2VsbHNcbiAgfSBlbHNlIHtcbiAgICBjZWxscy5zb3J0KGNvbXBhcmVDZWxscylcbiAgICByZXR1cm4gY2VsbHNcbiAgfVxufVxuZXhwb3J0cy5ub3JtYWxpemUgPSBub3JtYWxpemVcblxuLy9SZW1vdmVzIGFsbCBkdXBsaWNhdGUgY2VsbHMgaW4gdGhlIGNvbXBsZXhcbmZ1bmN0aW9uIHVuaXF1ZShjZWxscykge1xuICBpZihjZWxscy5sZW5ndGggPT09IDApIHtcbiAgICByZXR1cm4gW11cbiAgfVxuICB2YXIgcHRyID0gMVxuICAgICwgbGVuID0gY2VsbHMubGVuZ3RoXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpKSB7XG4gICAgdmFyIGEgPSBjZWxsc1tpXVxuICAgIGlmKGNvbXBhcmVDZWxscyhhLCBjZWxsc1tpLTFdKSkge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBjZWxsc1twdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGNlbGxzLmxlbmd0aCA9IHB0clxuICByZXR1cm4gY2VsbHNcbn1cbmV4cG9ydHMudW5pcXVlID0gdW5pcXVlO1xuXG4vL0ZpbmRzIGEgY2VsbCBpbiBhIG5vcm1hbGl6ZWQgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBmaW5kQ2VsbChjZWxscywgYykge1xuICB2YXIgbG8gPSAwXG4gICAgLCBoaSA9IGNlbGxzLmxlbmd0aC0xXG4gICAgLCByICA9IC0xXG4gIHdoaWxlIChsbyA8PSBoaSkge1xuICAgIHZhciBtaWQgPSAobG8gKyBoaSkgPj4gMVxuICAgICAgLCBzICAgPSBjb21wYXJlQ2VsbHMoY2VsbHNbbWlkXSwgYylcbiAgICBpZihzIDw9IDApIHtcbiAgICAgIGlmKHMgPT09IDApIHtcbiAgICAgICAgciA9IG1pZFxuICAgICAgfVxuICAgICAgbG8gPSBtaWQgKyAxXG4gICAgfSBlbHNlIGlmKHMgPiAwKSB7XG4gICAgICBoaSA9IG1pZCAtIDFcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJcbn1cbmV4cG9ydHMuZmluZENlbGwgPSBmaW5kQ2VsbDtcblxuLy9CdWlsZHMgYW4gaW5kZXggZm9yIGFuIG4tY2VsbC4gIFRoaXMgaXMgbW9yZSBnZW5lcmFsIHRoYW4gZHVhbCwgYnV0IGxlc3MgZWZmaWNpZW50XG5mdW5jdGlvbiBpbmNpZGVuY2UoZnJvbV9jZWxscywgdG9fY2VsbHMpIHtcbiAgdmFyIGluZGV4ID0gbmV3IEFycmF5KGZyb21fY2VsbHMubGVuZ3RoKVxuICBmb3IodmFyIGk9MCwgaWw9aW5kZXgubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICBpbmRleFtpXSA9IFtdXG4gIH1cbiAgdmFyIGIgPSBbXVxuICBmb3IodmFyIGk9MCwgbj10b19jZWxscy5sZW5ndGg7IGk8bjsgKytpKSB7XG4gICAgdmFyIGMgPSB0b19jZWxsc1tpXVxuICAgIHZhciBjbCA9IGMubGVuZ3RoXG4gICAgZm9yKHZhciBrPTEsIGtuPSgxPDxjbCk7IGs8a247ICsraykge1xuICAgICAgYi5sZW5ndGggPSBiaXRzLnBvcENvdW50KGspXG4gICAgICB2YXIgbCA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPGNsOyArK2opIHtcbiAgICAgICAgaWYoayAmICgxPDxqKSkge1xuICAgICAgICAgIGJbbCsrXSA9IGNbal1cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgdmFyIGlkeD1maW5kQ2VsbChmcm9tX2NlbGxzLCBiKVxuICAgICAgaWYoaWR4IDwgMCkge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgd2hpbGUodHJ1ZSkge1xuICAgICAgICBpbmRleFtpZHgrK10ucHVzaChpKVxuICAgICAgICBpZihpZHggPj0gZnJvbV9jZWxscy5sZW5ndGggfHwgY29tcGFyZUNlbGxzKGZyb21fY2VsbHNbaWR4XSwgYikgIT09IDApIHtcbiAgICAgICAgICBicmVha1xuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICB9XG4gIHJldHVybiBpbmRleFxufVxuZXhwb3J0cy5pbmNpZGVuY2UgPSBpbmNpZGVuY2VcblxuLy9Db21wdXRlcyB0aGUgZHVhbCBvZiB0aGUgbWVzaC4gIFRoaXMgaXMgYmFzaWNhbGx5IGFuIG9wdGltaXplZCB2ZXJzaW9uIG9mIGJ1aWxkSW5kZXggZm9yIHRoZSBzaXR1YXRpb24gd2hlcmUgZnJvbV9jZWxscyBpcyBqdXN0IHRoZSBsaXN0IG9mIHZlcnRpY2VzXG5mdW5jdGlvbiBkdWFsKGNlbGxzLCB2ZXJ0ZXhfY291bnQpIHtcbiAgaWYoIXZlcnRleF9jb3VudCkge1xuICAgIHJldHVybiBpbmNpZGVuY2UodW5pcXVlKHNrZWxldG9uKGNlbGxzLCAwKSksIGNlbGxzLCAwKVxuICB9XG4gIHZhciByZXMgPSBuZXcgQXJyYXkodmVydGV4X2NvdW50KVxuICBmb3IodmFyIGk9MDsgaTx2ZXJ0ZXhfY291bnQ7ICsraSkge1xuICAgIHJlc1tpXSA9IFtdXG4gIH1cbiAgZm9yKHZhciBpPTAsIGxlbj1jZWxscy5sZW5ndGg7IGk8bGVuOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTAsIGNsPWMubGVuZ3RoOyBqPGNsOyArK2opIHtcbiAgICAgIHJlc1tjW2pdXS5wdXNoKGkpXG4gICAgfVxuICB9XG4gIHJldHVybiByZXNcbn1cbmV4cG9ydHMuZHVhbCA9IGR1YWxcblxuLy9FbnVtZXJhdGVzIGFsbCBjZWxscyBpbiB0aGUgY29tcGxleFxuZnVuY3Rpb24gZXhwbG9kZShjZWxscykge1xuICB2YXIgcmVzdWx0ID0gW11cbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgICAgLCBjbCA9IGMubGVuZ3RofDBcbiAgICBmb3IodmFyIGo9MSwgamw9KDE8PGNsKTsgajxqbDsgKytqKSB7XG4gICAgICB2YXIgYiA9IFtdXG4gICAgICBmb3IodmFyIGs9MDsgazxjbDsgKytrKSB7XG4gICAgICAgIGlmKChqID4+PiBrKSAmIDEpIHtcbiAgICAgICAgICBiLnB1c2goY1trXSlcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmVzdWx0LnB1c2goYilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIG5vcm1hbGl6ZShyZXN1bHQpXG59XG5leHBvcnRzLmV4cGxvZGUgPSBleHBsb2RlXG5cbi8vRW51bWVyYXRlcyBhbGwgb2YgdGhlIG4tY2VsbHMgb2YgYSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIHNrZWxldG9uKGNlbGxzLCBuKSB7XG4gIGlmKG4gPCAwKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cbiAgdmFyIHJlc3VsdCA9IFtdXG4gICAgLCBrMCAgICAgPSAoMTw8KG4rMSkpLTFcbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBrPWswOyBrPCgxPDxjLmxlbmd0aCk7IGs9Yml0cy5uZXh0Q29tYmluYXRpb24oaykpIHtcbiAgICAgIHZhciBiID0gbmV3IEFycmF5KG4rMSlcbiAgICAgICAgLCBsID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8Yy5sZW5ndGg7ICsraikge1xuICAgICAgICBpZihrICYgKDE8PGopKSB7XG4gICAgICAgICAgYltsKytdID0gY1tqXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXN1bHQucHVzaChiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gbm9ybWFsaXplKHJlc3VsdClcbn1cbmV4cG9ydHMuc2tlbGV0b24gPSBza2VsZXRvbjtcblxuLy9Db21wdXRlcyB0aGUgYm91bmRhcnkgb2YgYWxsIGNlbGxzLCBkb2VzIG5vdCByZW1vdmUgZHVwbGljYXRlc1xuZnVuY3Rpb24gYm91bmRhcnkoY2VsbHMpIHtcbiAgdmFyIHJlcyA9IFtdXG4gIGZvcih2YXIgaT0wLGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wLGNsPWMubGVuZ3RoOyBqPGNsOyArK2opIHtcbiAgICAgIHZhciBiID0gbmV3IEFycmF5KGMubGVuZ3RoLTEpXG4gICAgICBmb3IodmFyIGs9MCwgbD0wOyBrPGNsOyArK2spIHtcbiAgICAgICAgaWYoayAhPT0gaikge1xuICAgICAgICAgIGJbbCsrXSA9IGNba11cbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgcmVzLnB1c2goYilcbiAgICB9XG4gIH1cbiAgcmV0dXJuIG5vcm1hbGl6ZShyZXMpXG59XG5leHBvcnRzLmJvdW5kYXJ5ID0gYm91bmRhcnk7XG5cbi8vQ29tcHV0ZXMgY29ubmVjdGVkIGNvbXBvbmVudHMgZm9yIGEgZGVuc2UgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBjb25uZWN0ZWRDb21wb25lbnRzX2RlbnNlKGNlbGxzLCB2ZXJ0ZXhfY291bnQpIHtcbiAgdmFyIGxhYmVscyA9IG5ldyBVbmlvbkZpbmQodmVydGV4X2NvdW50KVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MDsgajxjLmxlbmd0aDsgKytqKSB7XG4gICAgICBmb3IodmFyIGs9aisxOyBrPGMubGVuZ3RoOyArK2spIHtcbiAgICAgICAgbGFiZWxzLmxpbmsoY1tqXSwgY1trXSlcbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgdmFyIGNvbXBvbmVudHMgPSBbXVxuICAgICwgY29tcG9uZW50X2xhYmVscyA9IGxhYmVscy5yYW5rc1xuICBmb3IodmFyIGk9MDsgaTxjb21wb25lbnRfbGFiZWxzLmxlbmd0aDsgKytpKSB7XG4gICAgY29tcG9uZW50X2xhYmVsc1tpXSA9IC0xXG4gIH1cbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgbCA9IGxhYmVscy5maW5kKGNlbGxzW2ldWzBdKVxuICAgIGlmKGNvbXBvbmVudF9sYWJlbHNbbF0gPCAwKSB7XG4gICAgICBjb21wb25lbnRfbGFiZWxzW2xdID0gY29tcG9uZW50cy5sZW5ndGhcbiAgICAgIGNvbXBvbmVudHMucHVzaChbY2VsbHNbaV0uc2xpY2UoMCldKVxuICAgIH0gZWxzZSB7XG4gICAgICBjb21wb25lbnRzW2NvbXBvbmVudF9sYWJlbHNbbF1dLnB1c2goY2VsbHNbaV0uc2xpY2UoMCkpXG4gICAgfVxuICB9XG4gIHJldHVybiBjb21wb25lbnRzXG59XG5cbi8vQ29tcHV0ZXMgY29ubmVjdGVkIGNvbXBvbmVudHMgZm9yIGEgc3BhcnNlIGdyYXBoXG5mdW5jdGlvbiBjb25uZWN0ZWRDb21wb25lbnRzX3NwYXJzZShjZWxscykge1xuICB2YXIgdmVydGljZXMgID0gdW5pcXVlKG5vcm1hbGl6ZShza2VsZXRvbihjZWxscywgMCkpKVxuICAgICwgbGFiZWxzICAgID0gbmV3IFVuaW9uRmluZCh2ZXJ0aWNlcy5sZW5ndGgpXG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wOyBqPGMubGVuZ3RoOyArK2opIHtcbiAgICAgIHZhciB2aiA9IGZpbmRDZWxsKHZlcnRpY2VzLCBbY1tqXV0pXG4gICAgICBmb3IodmFyIGs9aisxOyBrPGMubGVuZ3RoOyArK2spIHtcbiAgICAgICAgbGFiZWxzLmxpbmsodmosIGZpbmRDZWxsKHZlcnRpY2VzLCBbY1trXV0pKVxuICAgICAgfVxuICAgIH1cbiAgfVxuICB2YXIgY29tcG9uZW50cyAgICAgICAgPSBbXVxuICAgICwgY29tcG9uZW50X2xhYmVscyAgPSBsYWJlbHMucmFua3NcbiAgZm9yKHZhciBpPTA7IGk8Y29tcG9uZW50X2xhYmVscy5sZW5ndGg7ICsraSkge1xuICAgIGNvbXBvbmVudF9sYWJlbHNbaV0gPSAtMVxuICB9XG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGwgPSBsYWJlbHMuZmluZChmaW5kQ2VsbCh2ZXJ0aWNlcywgW2NlbGxzW2ldWzBdXSkpO1xuICAgIGlmKGNvbXBvbmVudF9sYWJlbHNbbF0gPCAwKSB7XG4gICAgICBjb21wb25lbnRfbGFiZWxzW2xdID0gY29tcG9uZW50cy5sZW5ndGhcbiAgICAgIGNvbXBvbmVudHMucHVzaChbY2VsbHNbaV0uc2xpY2UoMCldKVxuICAgIH0gZWxzZSB7XG4gICAgICBjb21wb25lbnRzW2NvbXBvbmVudF9sYWJlbHNbbF1dLnB1c2goY2VsbHNbaV0uc2xpY2UoMCkpXG4gICAgfVxuICB9XG4gIHJldHVybiBjb21wb25lbnRzXG59XG5cbi8vQ29tcHV0ZXMgY29ubmVjdGVkIGNvbXBvbmVudHMgZm9yIGEgY2VsbCBjb21wbGV4XG5mdW5jdGlvbiBjb25uZWN0ZWRDb21wb25lbnRzKGNlbGxzLCB2ZXJ0ZXhfY291bnQpIHtcbiAgaWYodmVydGV4X2NvdW50KSB7XG4gICAgcmV0dXJuIGNvbm5lY3RlZENvbXBvbmVudHNfZGVuc2UoY2VsbHMsIHZlcnRleF9jb3VudClcbiAgfVxuICByZXR1cm4gY29ubmVjdGVkQ29tcG9uZW50c19zcGFyc2UoY2VsbHMpXG59XG5leHBvcnRzLmNvbm5lY3RlZENvbXBvbmVudHMgPSBjb25uZWN0ZWRDb21wb25lbnRzXG4iLCJcInVzZSBzdHJpY3RcIlxuXG5mdW5jdGlvbiB1bmlxdWVfcHJlZChsaXN0LCBjb21wYXJlKSB7XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBsaXN0Lmxlbmd0aFxuICAgICwgYT1saXN0WzBdLCBiPWxpc3RbMF1cbiAgZm9yKHZhciBpPTE7IGk8bGVuOyArK2kpIHtcbiAgICBiID0gYVxuICAgIGEgPSBsaXN0W2ldXG4gICAgaWYoY29tcGFyZShhLCBiKSkge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBsaXN0W3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgbGlzdC5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGxpc3Rcbn1cblxuZnVuY3Rpb24gdW5pcXVlX2VxKGxpc3QpIHtcbiAgdmFyIHB0ciA9IDFcbiAgICAsIGxlbiA9IGxpc3QubGVuZ3RoXG4gICAgLCBhPWxpc3RbMF0sIGIgPSBsaXN0WzBdXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpLCBiPWEpIHtcbiAgICBiID0gYVxuICAgIGEgPSBsaXN0W2ldXG4gICAgaWYoYSAhPT0gYikge1xuICAgICAgaWYoaSA9PT0gcHRyKSB7XG4gICAgICAgIHB0cisrXG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICBsaXN0W3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgbGlzdC5sZW5ndGggPSBwdHJcbiAgcmV0dXJuIGxpc3Rcbn1cblxuZnVuY3Rpb24gdW5pcXVlKGxpc3QsIGNvbXBhcmUsIHNvcnRlZCkge1xuICBpZihsaXN0Lmxlbmd0aCA9PT0gMCkge1xuICAgIHJldHVybiBsaXN0XG4gIH1cbiAgaWYoY29tcGFyZSkge1xuICAgIGlmKCFzb3J0ZWQpIHtcbiAgICAgIGxpc3Quc29ydChjb21wYXJlKVxuICAgIH1cbiAgICByZXR1cm4gdW5pcXVlX3ByZWQobGlzdCwgY29tcGFyZSlcbiAgfVxuICBpZighc29ydGVkKSB7XG4gICAgbGlzdC5zb3J0KClcbiAgfVxuICByZXR1cm4gdW5pcXVlX2VxKGxpc3QpXG59XG5cbm1vZHVsZS5leHBvcnRzID0gdW5pcXVlXG4iLCJcInVzZSBzdHJpY3RcIlxuXG52YXIgY2ggPSByZXF1aXJlKFwiaW5jcmVtZW50YWwtY29udmV4LWh1bGxcIilcbnZhciB1bmlxID0gcmVxdWlyZShcInVuaXFcIilcblxubW9kdWxlLmV4cG9ydHMgPSB0cmlhbmd1bGF0ZVxuXG5mdW5jdGlvbiBMaWZ0ZWRQb2ludChwLCBpKSB7XG4gIHRoaXMucG9pbnQgPSBwXG4gIHRoaXMuaW5kZXggPSBpXG59XG5cbmZ1bmN0aW9uIGNvbXBhcmVMaWZ0ZWQoYSwgYikge1xuICB2YXIgYXAgPSBhLnBvaW50XG4gIHZhciBicCA9IGIucG9pbnRcbiAgdmFyIGQgPSBhcC5sZW5ndGhcbiAgZm9yKHZhciBpPTA7IGk8ZDsgKytpKSB7XG4gICAgdmFyIHMgPSBicFtpXSAtIGFwW2ldXG4gICAgaWYocykge1xuICAgICAgcmV0dXJuIHNcbiAgICB9XG4gIH1cbiAgcmV0dXJuIDBcbn1cblxuZnVuY3Rpb24gdHJpYW5ndWxhdGUxRChuLCBwb2ludHMsIGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgaWYobiA9PT0gMSkge1xuICAgIGlmKGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgICAgIHJldHVybiBbIFstMSwgMF0gXVxuICAgIH0gZWxzZSB7XG4gICAgICByZXR1cm4gW11cbiAgICB9XG4gIH1cbiAgdmFyIGxpZnRlZCA9IHBvaW50cy5tYXAoZnVuY3Rpb24ocCwgaSkge1xuICAgIHJldHVybiBbIHBbMF0sIGkgXVxuICB9KVxuICBsaWZ0ZWQuc29ydChmdW5jdGlvbihhLGIpIHtcbiAgICByZXR1cm4gYVswXSAtIGJbMF1cbiAgfSlcbiAgdmFyIGNlbGxzID0gbmV3IEFycmF5KG4gLSAxKVxuICBmb3IodmFyIGk9MTsgaTxuOyArK2kpIHtcbiAgICB2YXIgYSA9IGxpZnRlZFtpLTFdXG4gICAgdmFyIGIgPSBsaWZ0ZWRbaV1cbiAgICBjZWxsc1tpLTFdID0gWyBhWzFdLCBiWzFdIF1cbiAgfVxuICBpZihpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gICAgY2VsbHMucHVzaChcbiAgICAgIFsgLTEsIGNlbGxzWzBdWzFdLCBdLFxuICAgICAgWyBjZWxsc1tuLTFdWzFdLCAtMSBdKVxuICB9XG4gIHJldHVybiBjZWxsc1xufVxuXG5mdW5jdGlvbiB0cmlhbmd1bGF0ZShwb2ludHMsIGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgdmFyIG4gPSBwb2ludHMubGVuZ3RoXG4gIGlmKG4gPT09IDApIHtcbiAgICByZXR1cm4gW11cbiAgfVxuICBcbiAgdmFyIGQgPSBwb2ludHNbMF0ubGVuZ3RoXG4gIGlmKGQgPCAxKSB7XG4gICAgcmV0dXJuIFtdXG4gIH1cblxuICAvL1NwZWNpYWwgY2FzZTogIEZvciAxRCB3ZSBjYW4ganVzdCBzb3J0IHRoZSBwb2ludHNcbiAgaWYoZCA9PT0gMSkge1xuICAgIHJldHVybiB0cmlhbmd1bGF0ZTFEKG4sIHBvaW50cywgaW5jbHVkZVBvaW50QXRJbmZpbml0eSlcbiAgfVxuICBcbiAgLy9MaWZ0IHBvaW50cywgc29ydFxuICB2YXIgbGlmdGVkID0gbmV3IEFycmF5KG4pXG4gIHZhciB1cHBlciA9IDEuMFxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICB2YXIgcCA9IHBvaW50c1tpXVxuICAgIHZhciB4ID0gbmV3IEFycmF5KGQrMSlcbiAgICB2YXIgbCA9IDAuMFxuICAgIGZvcih2YXIgaj0wOyBqPGQ7ICsraikge1xuICAgICAgdmFyIHYgPSBwW2pdXG4gICAgICB4W2pdID0gdlxuICAgICAgbCArPSB2ICogdlxuICAgIH1cbiAgICB4W2RdID0gbFxuICAgIGxpZnRlZFtpXSA9IG5ldyBMaWZ0ZWRQb2ludCh4LCBpKVxuICAgIHVwcGVyID0gTWF0aC5tYXgobCwgdXBwZXIpXG4gIH1cbiAgdW5pcShsaWZ0ZWQsIGNvbXBhcmVMaWZ0ZWQpXG4gIFxuICAvL0RvdWJsZSBwb2ludHNcbiAgbiA9IGxpZnRlZC5sZW5ndGhcblxuICAvL0NyZWF0ZSBuZXcgbGlzdCBvZiBwb2ludHNcbiAgdmFyIGRwb2ludHMgPSBuZXcgQXJyYXkobiArIGQgKyAxKVxuICB2YXIgZGluZGV4ID0gbmV3IEFycmF5KG4gKyBkICsgMSlcblxuICAvL0FkZCBzdGVpbmVyIHBvaW50cyBhdCB0b3BcbiAgdmFyIHUgPSAoZCsxKSAqIChkKzEpICogdXBwZXJcbiAgdmFyIHkgPSBuZXcgQXJyYXkoZCsxKVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgeVtpXSA9IDAuMFxuICB9XG4gIHlbZF0gPSB1XG5cbiAgZHBvaW50c1swXSA9IHkuc2xpY2UoKVxuICBkaW5kZXhbMF0gPSAtMVxuXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB2YXIgeCA9IHkuc2xpY2UoKVxuICAgIHhbaV0gPSAxXG4gICAgZHBvaW50c1tpKzFdID0geFxuICAgIGRpbmRleFtpKzFdID0gLTFcbiAgfVxuXG4gIC8vQ29weSByZXN0IG9mIHRoZSBwb2ludHMgb3ZlclxuICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICB2YXIgaCA9IGxpZnRlZFtpXVxuICAgIGRwb2ludHNbaSArIGQgKyAxXSA9IGgucG9pbnRcbiAgICBkaW5kZXhbaSArIGQgKyAxXSA9ICBoLmluZGV4XG4gIH1cblxuICAvL0NvbnN0cnVjdCBjb252ZXggaHVsbFxuICB2YXIgaHVsbCA9IGNoKGRwb2ludHMsIGZhbHNlKVxuICBpZihpbmNsdWRlUG9pbnRBdEluZmluaXR5KSB7XG4gICAgaHVsbCA9IGh1bGwuZmlsdGVyKGZ1bmN0aW9uKGNlbGwpIHtcbiAgICAgIHZhciBjb3VudCA9IDBcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgdmFyIHYgPSBkaW5kZXhbY2VsbFtqXV1cbiAgICAgICAgaWYodiA8IDApIHtcbiAgICAgICAgICBpZigrK2NvdW50ID49IDIpIHtcbiAgICAgICAgICAgIHJldHVybiBmYWxzZVxuICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBjZWxsW2pdID0gdlxuICAgICAgfVxuICAgICAgcmV0dXJuIHRydWVcbiAgICB9KVxuICB9IGVsc2Uge1xuICAgIGh1bGwgPSBodWxsLmZpbHRlcihmdW5jdGlvbihjZWxsKSB7XG4gICAgICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgICAgIHZhciB2ID0gZGluZGV4W2NlbGxbaV1dXG4gICAgICAgIGlmKHYgPCAwKSB7XG4gICAgICAgICAgcmV0dXJuIGZhbHNlXG4gICAgICAgIH1cbiAgICAgICAgY2VsbFtpXSA9IHZcbiAgICAgIH1cbiAgICAgIHJldHVybiB0cnVlXG4gICAgfSlcbiAgfVxuXG4gIGlmKGQgJiAxKSB7XG4gICAgZm9yKHZhciBpPTA7IGk8aHVsbC5sZW5ndGg7ICsraSkge1xuICAgICAgdmFyIGggPSBodWxsW2ldXG4gICAgICB2YXIgeCA9IGhbMF1cbiAgICAgIGhbMF0gPSBoWzFdXG4gICAgICBoWzFdID0geFxuICAgIH1cbiAgfVxuXG4gIHJldHVybiBodWxsXG59IiwiXG5tb2R1bGUuZXhwb3J0cyA9IHBhcnNlXG5cbi8qKlxuICogZXhwZWN0ZWQgYXJndW1lbnQgbGVuZ3Roc1xuICogQHR5cGUge09iamVjdH1cbiAqL1xuXG52YXIgbGVuZ3RoID0ge2E6IDcsIGM6IDYsIGg6IDEsIGw6IDIsIG06IDIsIHE6IDQsIHM6IDQsIHQ6IDIsIHY6IDEsIHo6IDB9XG5cbi8qKlxuICogc2VnbWVudCBwYXR0ZXJuXG4gKiBAdHlwZSB7UmVnRXhwfVxuICovXG5cbnZhciBzZWdtZW50ID0gLyhbYXN0dnpxbWhsY10pKFteYXN0dnpxbWhsY10qKS9pZ1xuXG4vKipcbiAqIHBhcnNlIGFuIHN2ZyBwYXRoIGRhdGEgc3RyaW5nLiBHZW5lcmF0ZXMgYW4gQXJyYXlcbiAqIG9mIGNvbW1hbmRzIHdoZXJlIGVhY2ggY29tbWFuZCBpcyBhbiBBcnJheSBvZiB0aGVcbiAqIGZvcm0gYFtjb21tYW5kLCBhcmcxLCBhcmcyLCAuLi5dYFxuICpcbiAqIEBwYXJhbSB7U3RyaW5nfSBwYXRoXG4gKiBAcmV0dXJuIHtBcnJheX1cbiAqL1xuXG5mdW5jdGlvbiBwYXJzZShwYXRoKSB7XG5cdHZhciBkYXRhID0gW11cblx0cGF0aC5yZXBsYWNlKHNlZ21lbnQsIGZ1bmN0aW9uKF8sIGNvbW1hbmQsIGFyZ3Mpe1xuXHRcdHZhciB0eXBlID0gY29tbWFuZC50b0xvd2VyQ2FzZSgpXG5cdFx0YXJncyA9IHBhcnNlVmFsdWVzKGFyZ3MpXG5cblx0XHQvLyBvdmVybG9hZGVkIG1vdmVUb1xuXHRcdGlmICh0eXBlID09ICdtJyAmJiBhcmdzLmxlbmd0aCA+IDIpIHtcblx0XHRcdGRhdGEucHVzaChbY29tbWFuZF0uY29uY2F0KGFyZ3Muc3BsaWNlKDAsIDIpKSlcblx0XHRcdHR5cGUgPSAnbCdcblx0XHRcdGNvbW1hbmQgPSBjb21tYW5kID09ICdtJyA/ICdsJyA6ICdMJ1xuXHRcdH1cblxuXHRcdHdoaWxlICh0cnVlKSB7XG5cdFx0XHRpZiAoYXJncy5sZW5ndGggPT0gbGVuZ3RoW3R5cGVdKSB7XG5cdFx0XHRcdGFyZ3MudW5zaGlmdChjb21tYW5kKVxuXHRcdFx0XHRyZXR1cm4gZGF0YS5wdXNoKGFyZ3MpXG5cdFx0XHR9XG5cdFx0XHRpZiAoYXJncy5sZW5ndGggPCBsZW5ndGhbdHlwZV0pIHRocm93IG5ldyBFcnJvcignbWFsZm9ybWVkIHBhdGggZGF0YScpXG5cdFx0XHRkYXRhLnB1c2goW2NvbW1hbmRdLmNvbmNhdChhcmdzLnNwbGljZSgwLCBsZW5ndGhbdHlwZV0pKSlcblx0XHR9XG5cdH0pXG5cdHJldHVybiBkYXRhXG59XG5cbmZ1bmN0aW9uIHBhcnNlVmFsdWVzKGFyZ3Mpe1xuXHRhcmdzID0gYXJncy5tYXRjaCgvLT9bLjAtOV0rKD86ZVstK10/XFxkKyk/L2lnKVxuXHRyZXR1cm4gYXJncyA/IGFyZ3MubWFwKE51bWJlcikgOiBbXVxufVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgc2lnbiA9IHJlcXVpcmUoJy4vdXRpbGl0aWVzLmpzJykuc2lnbjtcbnZhciBjYWxjdWxhdGVEaXN0YW5jZSA9IHJlcXVpcmUoJy4vdXRpbGl0aWVzLmpzJykuZGlzdGFuY2U7XG5cbi8vIHZhciBwb2ludHMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5wb2ludHM7XG4vLyB2YXIgY2l0eVNldCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLmNpdHlTZXQ7XG4vLyB2YXIgdGV4dFBvaW50c0lkID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykudGV4dFBvaW50c0lkO1xuLy8gdmFyIHBvc3NpYmxlU3RhcnRQb2ludHNJZCA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLnBvc3NpYmxlU3RhcnRQb2ludHNJZDtcblxudmFyIGxpdmVNb3VzZVBvc2l0aW9uID0gcmVxdWlyZSgnLi9tb3VzZS5qcycpO1xuXG52YXIgVmVjdG9yID0gcmVxdWlyZSgnLi92ZWN0b3IuanMnKTtcblxudmFyIHJhbmRvbSA9IE1hdGgucmFuZG9tO1xudmFyIGZsb29yID0gTWF0aC5mbG9vcjtcbnZhciBXRUlHSFQgPSAxMDtcbi8vIHZhciBSRVBVTFNJT04gPSAwLjA1O1xuLy8gdmFyIFJFUFVMU0lPTlNQRUVEID0gMC4wMDI7XG4vLyB2YXIgQU5UVkVMT0NJVFkgPSAwLjAwMTtcblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihjb250YWluZXIsIGluaXRQb2ludHMsIG9wdGlvbnMpe1xuXG4gICAgY29uc29sZS5sb2coJ09wdGlvbnMgYW50IDonLCBvcHRpb25zKTtcbiAgICAvLyBEZWZpbmUgdGhvc2UgcGFyYW1ldGVycyBhcyBhdHRyaWJ1dGVzIG9mIEFudCBvYmplY3QgP1xuICAgIHZhciBSRVBVTFNJT04gPSBvcHRpb25zLnJlcFNpemU7XG4gICAgdmFyIFJFUFVMU0lPTlNQRUVEID0gb3B0aW9ucy5yZXBTcGVlZDtcbiAgICB2YXIgQU5UVkVMT0NJVFkgPSBvcHRpb25zLnZlbG9jaXR5O1xuICAgIHZhciBJTlRFTExJR0VOQ0UgPSBvcHRpb25zLmludGVsbGlnZW5jZTtcblxuICAgIHZhciBtb3VzZSA9IGxpdmVNb3VzZVBvc2l0aW9uKGNvbnRhaW5lcik7XG5cbiAgICB2YXIgcG9pbnRzID0gaW5pdFBvaW50cy5wb2ludHM7XG4gICAgdmFyIGNpdHlTZXQgPSBpbml0UG9pbnRzLmNpdHlTZXQ7XG4gICAgdmFyIHRleHRQb2ludHNJZCA9IGluaXRQb2ludHMudGV4dFBvaW50c0lkO1xuICAgIHZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSBpbml0UG9pbnRzLnBvc3NpYmxlU3RhcnRQb2ludHNJZDtcblxuXG4gICAgZnVuY3Rpb24gQW50KHBvaW50KSB7XG4gICAgICAgIHRoaXMueCA9IHBvaW50Lng7ICAgICAgICAgICAgICAgIFxuICAgICAgICB0aGlzLnkgPSBwb2ludC55O1xuICAgICAgICB0aGlzLnZlbG9jaXR5ID0gQU5UVkVMT0NJVFk7XG4gICAgICAgIHRoaXMuaW50ZWxsaWdlbmNlID0gSU5URUxMSUdFTkNFO1xuICAgICAgICB0aGlzLnJlcFNpemUgPSBSRVBVTFNJT047XG4gICAgICAgIHRoaXMucmVwU3BlZWQgPSBSRVBVTFNJT05TUEVFRDtcbiAgICAgICAgdGhpcy5lZGdlID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLnN0YXRlID0gXCJmb3JhZ2VcIjtcbiAgICAgICAgdGhpcy5lZGdlcyA9IFtdO1xuICAgICAgICB0aGlzLmxhc3RDaXR5ID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLm9yaWdpbiA9IHBvaW50O1xuICAgICAgICB0aGlzLmRlc3RpbmF0aW9uID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLm9yaWVudGF0aW9uID0gdW5kZWZpbmVkO1xuICAgICAgICB0aGlzLmRpcmVjdGlvbiA9IG5ldyBWZWN0b3IoMCwwKTtcbiAgICAgICAgdGhpcy5wcm9nID0gMDtcbiAgICB9XG4gICAgLy8gZm9yYWdlOiB0aGUgYW50IHdhbmRlcnMgYXJvdW5kIHdpdGhvdXQgYW55IHBoZXJvbW9uIGRlcG9zaXRpb25cbiAgICAvLyBvbmNlIGl0IGZpbmRzIGEgY2l0eSwgaXQgc3RhcnRzIHJlbWVtYmVyaW5nIHRoZSBub2RlcyBpdCBnb2VzIHRocm91Z2hcbiAgICAvLyB3aGVuIGl0IGZpbmRzIGFub3RoZXIgY2l0eSwgaXQgY29tcHV0ZXMgdGhlIHBhdGggbGVuZ3RoIGFuZCBhZGRzIHBoZXJvbW9ucyBvbmUgZWFjaCBlZGdlc1xuICAgIC8vIHByb3BvcnRpb25uYWx5IHRvIHRoZSBzaG9ydGVzdG5lc3Mgb2YgdGhlIHBhdGhcbiAgICAvLyBpdCByZXNldHMgdGhlIGxpc3Qgb2Ygbm9kZXMgYW5kIGNvbnRpbnVlc1xuICAgIC8vIHdoaWxlIGZvcmFnaW5nIHRoZSBhbnQgY2hvc2VzIHRoZSBwYXRoIHdpdGggYSBwaGVyb21vbiBwcmVmZXJlbmNlXG5cblxuICAgIC8vIHN0YXRpYyBtZXRob2RzXG4gICAgQW50LmdlbmVyYXRlUmFuZFN0YXJ0UG9pbnQgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgdmFyIHJhbmRJZCA9IE1hdGguZmxvb3IocG9zc2libGVTdGFydFBvaW50c0lkLmxlbmd0aCAqIHJhbmRvbSgpKTtcbiAgICAgICAgdmFyIHJhbmRTdGFydFBvaW50ID0gcG9pbnRzW3Bvc3NpYmxlU3RhcnRQb2ludHNJZFtyYW5kSWRdXTtcbiAgICAgICAgcmV0dXJuIHJhbmRTdGFydFBvaW50O1xuICAgIH1cblxuXG4gICAgLy8gbWV0aG9kc1xuICAgIEFudC5wcm90b3R5cGUgPSB7XG5cbiAgICAgICAgdHJhbnNpdDogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHN3aXRjaCAodGhpcy5zdGF0ZSkge1xuICAgICAgICAgICAgY2FzZSBcImZvcmFnZVwiOlxuICAgICAgICAgICAgICAgIHZhciByZXMgPSB0aGlzLm1vdmUoKTtcbiAgICAgICAgICAgICAgICBpZiAocmVzLmNpdHlSZWFjaGVkKSB7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMuc3RhdGUgPSBcInBoZXJvbW9uaW5nXCI7XG4gICAgICAgICAgICAgICAgICAgIHRoaXMubGFzdENpdHkgPSB0aGlzLm9yaWdpbi5pZDtcbiAgICAgICAgICAgICAgICB9O1xuICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgY2FzZSBcInBoZXJvbW9uaW5nXCI6XG4gICAgICAgICAgICAgICAgdmFyIHJlcyA9IHRoaXMubW92ZSgpO1xuICAgICAgICAgICAgICAgIGlmIChyZXMuZWRnZUNoYW5nZWQpIHtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5lZGdlcy5wdXNoKHRoaXMuZWRnZSk7XG4gICAgICAgICAgICAgICAgICAgIC8vIGZvdW5kIGEgY2l0eVxuICAgICAgICAgICAgICAgICAgICBpZiAocmVzLmNpdHlSZWFjaGVkICYmICh0aGlzLm9yaWdpbi5pZCAhPSB0aGlzLmxhc3RDaXR5KSApe1xuICAgICAgICAgICAgICAgICAgICAgICAgLy8gY29tcHV0ZSB0aGUgbGVuZ3RoIG9mIHRoZSBwYXRoXG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgcGF0aExlbmd0aCA9IHRoaXMuZWRnZXMubWFwKGZ1bmN0aW9uKGUpe3JldHVybiBlLmRpc3RhbmNlfSkucmVkdWNlKGZ1bmN0aW9uKGEsYil7cmV0dXJuIGEgKyBifSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgZGVsdGFQaGVyb21vbmUgPSAxL3BhdGhMZW5ndGg7XG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzLmZvckVhY2goZnVuY3Rpb24oZSl7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGEgPSBlLnB0MSwgYiA9IGUucHQyLCB3ZWlnaHQgPSAxOyAgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgLy8gaW5jcmVhc2VkIGRyb3BwZWQgcGhlcm9tb25zIGZvciB0ZXh0RWRnZXNcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoKGNpdHlTZXQuaW5kZXhPZihhLmlkKSAhPSAtMSkgJiYgY2l0eVNldC5pbmRleE9mKGIuaWQpICE9IC0xICYmIChNYXRoLmFicyhhLmlkIC0gYi5pZCkgPT0gMSkpXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB3ZWlnaHQgKj0gV0VJR0hUO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBlLnBoZXJvbW9uICs9IChkZWx0YVBoZXJvbW9uZSAqIHdlaWdodCk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9KTtcblxuICAgICAgICAgICAgICAgICAgICAgICAgdGhpcy5lZGdlcyA9IFt0aGlzLmVkZ2VdO1xuICAgICAgICAgICAgICAgICAgICAgICAgdGhpcy5sYXN0Q2l0eSA9IHRoaXMub3JpZ2luLmlkO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgIH1cblxuICAgICAgICB9LFxuXG4gICAgICAgIHNldERpcmVjdGlvbjogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHZhciBwb3NzaWJsZUVkZ2VzID0gW107XG5cbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdGhpcy5vcmlnaW4ubmV4dHMubGVuZ3RoOyBpKyspXG4gICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgcG9zc2libGVFZGdlc1tpXSA9IHRoaXMub3JpZ2luLm5leHRzW2ldO1xuICAgICAgICAgICAgfSBcblxuICAgICAgICAgICAgLy8gY29uc29sZS5sb2coJ3NtZWxsczE6ICcsIHBvc3NpYmxlRWRnZXMpO1xuXG4gICAgICAgICAgICBwb3NzaWJsZUVkZ2VzLnNwbGljZShwb3NzaWJsZUVkZ2VzLmluZGV4T2YodGhpcy5lZGdlKSwxKTtcblxuICAgICAgICAgICAgLy8gZmxpcCBhIGNvaW4gYW5kIGVpdGhlciB0YWtlIHRoZSBzbWVsbGllc3QgcGF0aCBvciBhIHJhbmRvbSBvbmVcbiAgICAgICAgICAgIGlmIChyYW5kb20oKSA8IHRoaXMuaW50ZWxsaWdlbmNlKXtcbiAgICAgICAgICAgICAgICB2YXIgc21lbGxzID0gcG9zc2libGVFZGdlcy5tYXAoZnVuY3Rpb24oZSl7cmV0dXJuIGUucGhlcm9tb247fSk7XG4gICAgICAgICAgICAgICAgdmFyIGluZGV4ID0gc21lbGxzLmluZGV4T2YoTWF0aC5tYXguYXBwbHkoTWF0aCwgc21lbGxzKSk7XG4gICAgICAgICAgICAgICAgdGhpcy5lZGdlID0gcG9zc2libGVFZGdlc1tpbmRleF07XG4gICAgICAgICAgICB9IFxuICAgICAgICAgICAgZWxzZXtcbiAgICAgICAgICAgICAgICB0aGlzLmVkZ2UgPSBwb3NzaWJsZUVkZ2VzW2Zsb29yKHJhbmRvbSgpKnBvc3NpYmxlRWRnZXMubGVuZ3RoKV07XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgXG5cbiAgICAgICAgICAgIC8vIHNldCB0aGUgZGVzdGluYXRpb24gcG9pbnQsIGJlaW5nIGVkZ2UucHQxIG9yIGVkZ2UucHQyXG4gICAgICAgICAgICB0aGlzLmRlc3RpbmF0aW9uID0gKHRoaXMub3JpZ2luID09IHRoaXMuZWRnZS5wdDEpID8gdGhpcy5lZGdlLnB0MiA6IHRoaXMuZWRnZS5wdDE7XG5cbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLnggPSB0aGlzLmRlc3RpbmF0aW9uLnggLSB0aGlzLm9yaWdpbi54OyBcbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLnkgPSB0aGlzLmRlc3RpbmF0aW9uLnkgLSB0aGlzLm9yaWdpbi55O1xuXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi5ub3JtYWxpemUoKTtcbiAgICAgICAgfSxcblxuICAgICAgICBtb3ZlOiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgLy8gY29uc29sZS5sb2coJ21vdmUnKTtcbiAgICAgICAgICAgIHZhciBlZGdlQ2hhbmdlZDtcbiAgICAgICAgICAgIHZhciBjaXR5UmVhY2hlZCA9IGZhbHNlO1xuXG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi54ID0gdGhpcy5kZXN0aW5hdGlvbi54IC0gdGhpcy54OyBcbiAgICAgICAgICAgIHRoaXMuZGlyZWN0aW9uLnkgPSB0aGlzLmRlc3RpbmF0aW9uLnkgLSB0aGlzLnk7XG4gICAgICAgICAgICB0aGlzLmRpcmVjdGlvbi5ub3JtYWxpemUoKTtcblxuICAgICAgICAgICAgLy8gb24gZWRnZVxuICAgICAgICAgICAgaWYgKChjYWxjdWxhdGVEaXN0YW5jZSh0aGlzLCB0aGlzLmRlc3RpbmF0aW9uKSA+IHRoaXMucmVwU3BlZWQpKXtcblxuICAgICAgICAgICAgICAgIC8vIGEgZGVsdGEgbW92ZW1lbnQgd2lsbCBiZSBhcHBsaWVkIGlmIGNvbGxpc2lvbiB3aXRoIG9ic3RhY2xlIGRldGVjdGVkXG4gICAgICAgICAgICAgICAgdmFyIGRlbHRhID0gdGhpcy5hdm9pZE9ic3RhY2xlKCk7XG5cbiAgICAgICAgICAgICAgICB0aGlzLnggKz0gdGhpcy52ZWxvY2l0eSAqIHRoaXMuZGlyZWN0aW9uLnggKyBkZWx0YS54ICogdGhpcy5yZXBTcGVlZDtcbiAgICAgICAgICAgICAgICB0aGlzLnkgKz0gdGhpcy52ZWxvY2l0eSAqIHRoaXMuZGlyZWN0aW9uLnkgKyBkZWx0YS55ICogdGhpcy5yZXBTcGVlZDtcblxuICAgICAgICAgICAgICAgIHRoaXMucHJvZyA9IHRoaXMuY2FsY3VsYXRlUHJvZ3Jlc3Npb24oKTtcbiAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICBlZGdlQ2hhbmdlZCA9IGZhbHNlO1xuXG4gICAgICAgICAgICAvLyBvbiB2ZXJ0ZXhcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgLy8gY29uc29sZS5sb2coJ3JlYWNoZWQnKTtcbiAgICAgICAgICAgICAgICB0aGlzLnN0ZXAgPSAwO1xuICAgICAgICAgICAgICAgIHRoaXMucHJvZyA9IDA7XG4gICAgICAgICAgICAgICAgdGhpcy5vcmlnaW4gPSB0aGlzLmRlc3RpbmF0aW9uO1xuICAgICAgICAgICAgICAgIHRoaXMueCA9IHRoaXMub3JpZ2luLng7XG4gICAgICAgICAgICAgICAgdGhpcy55ID0gdGhpcy5vcmlnaW4ueTtcblxuICAgICAgICAgICAgICAgIHRoaXMuc2V0RGlyZWN0aW9uKCk7XG5cbiAgICAgICAgICAgICAgICBjaXR5UmVhY2hlZCA9IChjaXR5U2V0LmluZGV4T2YodGhpcy5vcmlnaW4uaWQpICE9IC0xKTtcbiAgICAgICAgICAgICAgICBlZGdlQ2hhbmdlZCA9IHRydWU7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXR1cm4ge2NpdHlSZWFjaGVkOiBjaXR5UmVhY2hlZCwgZWRnZUNoYW5nZWQ6IGVkZ2VDaGFuZ2VkfTtcbiAgICAgICAgfSxcblxuICAgICAgICBhdm9pZE9ic3RhY2xlOiBmdW5jdGlvbigpe1xuICAgICAgICAgICAgdmFyIGRpc3RhbmNlID0gY2FsY3VsYXRlRGlzdGFuY2UodGhpcywgbW91c2UpO1xuICAgICAgICBcbiAgICAgICAgICAgIGlmIChkaXN0YW5jZSA8PSB0aGlzLnJlcFNpemUpIHtcblxuICAgICAgICAgICAgICAgIHJldHVybiB7XG4gICAgICAgICAgICAgICAgICAgIC8vIGRlbHRhIG1vdmVtZW50IGlzIGNvbXBvc2VkIG9mIGEgcmVwdWxzaW9uIGRlbHRhIGFuZCBhIGNpcmN1bGFyIGRlbHRhIFxuICAgICAgICAgICAgICAgICAgICB4OiAodGhpcy54IC0gbW91c2UueCkvZGlzdGFuY2UgKyAodGhpcy55IC0gbW91c2UueSkvZGlzdGFuY2UgKiAxLFxuICAgICAgICAgICAgICAgICAgICB5OiAodGhpcy55IC0gbW91c2UueSkvZGlzdGFuY2UgLSAodGhpcy54IC0gbW91c2UueCkvZGlzdGFuY2UgKiAxXG4gICAgICAgICAgICAgICAgfTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGVsc2VcbiAgICAgICAgICAgICAgICByZXR1cm4ge3g6MCwgeTowfTtcbiAgICAgICAgfSxcblxuICAgICAgICBjYWxjdWxhdGVQcm9ncmVzc2lvbjogZnVuY3Rpb24oKXtcbiAgICAgICAgICAgIHZhciB2ID0gbmV3IFZlY3Rvcih0aGlzLnggLSB0aGlzLm9yaWdpbi54LCB0aGlzLnkgLSB0aGlzLm9yaWdpbi55KTtcbiAgICAgICAgICAgIHZhciBub3JtID0gdi5ub3JtKCk7XG5cbiAgICAgICAgICAgIHZhciB0aGV0YSA9ICh2LnggKiB0aGlzLmVkZ2UuZGlyZWN0aW9uLnggKyB2LnkgKiB0aGlzLmVkZ2UuZGlyZWN0aW9uLnkpIC8gbm9ybTtcbiAgICAgICAgICAgIHZhciBwcm9nID0gbm9ybSAqIE1hdGguYWJzKHRoZXRhKTtcbiAgICAgICAgICAgIC8vIHJldHVybnMgbGVuZ3RoIG9mIHByb2plY3Rpb24gb24gZWRnZVxuICAgICAgICAgICAgcmV0dXJuIHByb2c7XG4gICAgICAgIH1cblxuICAgIH07XG4gICAgcmV0dXJuIEFudDtcbn1cblxuIiwiJ3VzZSBzdHJpY3QnXG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24gKEFudCkge1xuXG5cdHZhciBuYkFudHNQZXJTdGVwID0gMTAwO1xuXG5cdGZ1bmN0aW9uIGNyZWF0ZUdyb3VwKHBvcHVsYXRpb24pe1xuXHRcdGZvciAodmFyIGkgPSAwOyBpIDwgbmJBbnRzUGVyU3RlcDsgaSsrKSB7XG5cdFx0XHR2YXIgbmV3QW50ID0gbmV3IEFudChBbnQuZ2VuZXJhdGVSYW5kU3RhcnRQb2ludCgpKTtcblx0XHRcdG5ld0FudC5zZXREaXJlY3Rpb24oKTtcblx0XHRcdHBvcHVsYXRpb24ucHVzaChuZXdBbnQpO1xuXHRcdH1cblxuLy8gXHRcdGNvbnNvbGUubG9nKCdDcmVhdGVkIEFudHMgR3JvdXA6IFxcXG4vLyAoKyAnICsgbmJBbnRzUGVyU3RlcCArICcpID0+ICcgKyBwb3B1bGF0aW9uLmxlbmd0aCk7XG5cblx0XHRyZXR1cm4gcG9wdWxhdGlvbjtcblx0fVxuXG5cdGZ1bmN0aW9uIHJlbW92ZUdyb3VwKHBvcHVsYXRpb24sIG5iRGVhZCl7XG5cdFx0cG9wdWxhdGlvbiA9IHBvcHVsYXRpb24uc2xpY2UoMCwgcG9wdWxhdGlvbi5sZW5ndGggLSBuYkRlYWQpO1xuXG4vLyBcdFx0Y29uc29sZS5sb2coJ1JlbW92ZWQgQW50cyBHcm91cDogXFxcbi8vICgtICcgKyBuYkFudHNQZXJTdGVwICsgJykgPT4gJyArIHBvcHVsYXRpb24ubGVuZ3RoKTtcblxuXHRcdHJldHVybiBwb3B1bGF0aW9uO1xuXG5cdH1cblxuXHRyZXR1cm4ge1xuXHRcdGNyZWF0ZTogY3JlYXRlR3JvdXAsXG5cdFx0cmVtb3ZlOiByZW1vdmVHcm91cFxuXHR9O1xuXG59XG5cdCIsIid1c2Ugc3RyaWN0J1xuXG52YXIgZHQgPSByZXF1aXJlKFwiZGVsYXVuYXktdHJpYW5ndWxhdGVcIik7XG5cbnZhciByYW5nZSA9IHJlcXVpcmUoJy4vdXRpbGl0aWVzLmpzJykucmFuZ2U7XG5cbi8vIHZhciBwb2ludHMgPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS5wb2ludHM7XG52YXIgdGV4dE1lc2ggPSByZXF1aXJlKCcuL2luaXRpYWxpemVQb2ludHMuanMnKS50ZXh0TWVzaDtcbnZhciBjaXR5U2V0ID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykuY2l0eVNldDtcbnZhciBuYlJhbmRvbVBvaW50cyA9IHJlcXVpcmUoJy4vaW5pdGlhbGl6ZVBvaW50cy5qcycpLm5iUmFuZG9tUG9pbnRzO1xudmFyIGZvcmNlZEVkZ2VzID0gcmVxdWlyZSgnLi9pbml0aWFsaXplUG9pbnRzLmpzJykuZm9yY2VkRWRnZXM7XG5cbnZhciBFZGdlID0gcmVxdWlyZSgnLi9lZGdlLmpzJyk7XG5cblxubW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbihwb2ludHMpe1xuICAgIC8vIHRyaWFuZ3VsYXRlXG4gICAgdmFyIGNlbGxzID0gZHQocG9pbnRzLm1hcChmdW5jdGlvbihwKXtcbiAgICAgICAgcmV0dXJuIFtwLngsIHAueV07XG4gICAgfSkpO1xuXG4gICAgdmFyIGVkZ2VzID0gW107XG4gICAgdmFyIHBlcm11dGF0aW9ucyA9IFtbMCwxXSwgWzAsMl0sIFsxLDJdXTtcblxuICAgIC8vIGZvcmNlIHRoZSBlZGdlcyBvZiB0aGUgdGV4dCB0byBiZSBlZGdlcyBvZiB0aGUgZ3JhcGhcbiAgICBpZiAodGV4dE1lc2gpIHtcbiAgICAgICAgcmFuZ2UoMCwgcG9pbnRzLmxlbmd0aCAtIG5iUmFuZG9tUG9pbnRzKS5mb3JFYWNoKGZ1bmN0aW9uKGlkKXtcbiAgICAgICAgICAgIHZhciBkaXJlY3RMaW5rID0gZm9yY2VkRWRnZXNbaWRdO1xuICAgICAgICAgICAgdmFyIHRleHRFZGdlID0gRWRnZS5jcmVhdGUocG9pbnRzW2lkXSwgcG9pbnRzW2RpcmVjdExpbmtdKTtcbiAgICAgICAgICAgIGVkZ2VzLnB1c2godGV4dEVkZ2UpO1xuICAgICAgICAgICAgcG9pbnRzW2lkXS5uZXh0cy5wdXNoKHRleHRFZGdlKTtcbiAgICAgICAgfSlcbiAgICB9XG5cblxuICAgIGNlbGxzLmZvckVhY2goZnVuY3Rpb24oY2VsbCl7XG4gICAgICAgXG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgMzsgKytpKXsgIC8vIGZvciBlYWNoIHBvaW50LmlkIGxpc3RlZCBpbiBjdXJyZW50IGNlbGxcbiAgICAgICAgICAgIHZhciBwdCA9IHBvaW50c1tjZWxsW2ldXTtcblxuICAgICAgICAgICAgZm9yICh2YXIgaiA9IDE7IGogPCAzOyArK2opeyBcblxuICAgICAgICAgICAgICAgIHZhciBwdGogPSBwb2ludHNbY2VsbFsoIGkgKyBqICkgJSAzXV07IC8vIHBpY2sgb25lIG9mIHRoZSBvdGhlciAyIHBvaW50cyBvZiB0aGUgY2VsbFxuICAgICAgICAgICAgICAgIHZhciBuZXdFZGdlID0gdW5kZWZpbmVkO1xuXG4gICAgICAgICAgICAgICAgLy8gaWYgcHQgYWxyZWFkeSBoYXMgbmV4dEVkZ2VzIC4uLlxuICAgICAgICAgICAgICAgIGlmIChwdC5uZXh0cy5sZW5ndGggIT0gMCkge1xuICAgICAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgLy8gLi4uIGdldCB0aGUgcG9pbnRzIGNvcnJlc3BvbmRpbmcgLi4uXG4gICAgICAgICAgICAgICAgICAgIHZhciB0ZW1wUG9pbnRzID0gcHQubmV4dHMubWFwKGZ1bmN0aW9uKGUpe1xuICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIFtlLnB0MSwgZS5wdDJdO1xuICAgICAgICAgICAgICAgICAgICB9KS5yZWR1Y2UoZnVuY3Rpb24oYSwgYil7XG4gICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGEuY29uY2F0KGIpO1xuICAgICAgICAgICAgICAgICAgICB9KTtcblxuICAgICAgICAgICAgICAgICAgICAvLyAuLi4gYW5kIGNoZWNrIGlmIHB0aiBhbHJlYWR5IGlzIHBhcnQgb2YgdGhlIGV4aXN0aW5nIG5leHRFZGdlcy4gSWYgbm90LCBhZGQgdGhlIGVkZ2UuXG4gICAgICAgICAgICAgICAgICAgIGlmICh0ZW1wUG9pbnRzLmluZGV4T2YocHRqKSA9PSAtMSl7XG4gICAgICAgICAgICAgICAgICAgICAgICBuZXdFZGdlID0gRWRnZS5jcmVhdGUocHQsIHB0aik7XG4gICAgICAgICAgICAgICAgICAgICAgICBlZGdlcy5wdXNoKG5ld0VkZ2UpO1xuICAgICAgICAgICAgICAgICAgICAgICAgcHQubmV4dHMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgbmV3RWRnZSA9IEVkZ2UuY3JlYXRlKHB0LCBwdGopO1xuICAgICAgICAgICAgICAgICAgICBlZGdlcy5wdXNoKG5ld0VkZ2UpO1xuICAgICAgICAgICAgICAgICAgICBwdC5uZXh0cy5wdXNoKG5ld0VkZ2UpO1xuICAgICAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgICAgIC8vIGFkZCBhbHNvIHRoZSBlZGdlIHRvIHRoZSBlZGdlJ3Mgb3RoZXIgcG9pbnQncyBuZXh0RWRnZXNcbiAgICAgICAgICAgICAgICBpZiAobmV3RWRnZSAhPSB1bmRlZmluZWQpe1xuICAgICAgICAgICAgICAgICAgICBwdGoubmV4dHMucHVzaChuZXdFZGdlKTtcbiAgICAgICAgICAgICAgICB9ICAgICAgICAgXG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIC8vIGFkZCB0aGUgdGV4dEVkZ2VzIHRvIG5leHRFZGdlcyBtYXBcbiAgICAgICAgICAgIGlmICh0ZXh0TWVzaCAmJiAoY2l0eVNldC5pbmRleE9mKHB0KSAhPSAtMSkpIHtcbiAgICAgICAgICAgICAgICB2YXIgdGV4dEVkZ2UgPSBFZGdlLmNyZWF0ZShwdCwgcG9pbnRzW3B0LmlkICsgMV0pO1xuICAgICAgICAgICAgICAgIGVkZ2VzLnB1c2godGV4dEVkZ2UpO1xuICAgICAgICAgICAgICAgIHB0Lm5leHRzLnB1c2godGV4dEVkZ2UpO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgIH1cbiAgICB9KTtcblxuICAgIHJldHVybiBlZGdlcztcbn07IiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgc3FydCA9IE1hdGguc3FydDtcbnZhciBwb3cgPSBNYXRoLnBvdztcbnZhciBhYnMgPSBNYXRoLmFicztcbnZhciBhdGFuID0gTWF0aC5hdGFuO1xuXG52YXIgVmVjdG9yID0gcmVxdWlyZSgnLi92ZWN0b3IuanMnKTtcblxuXG5mdW5jdGlvbiBFZGdlKHB0QSwgcHRCKSB7XG4gICAgdmFyIGRpc3RhbmNlID0gc3FydCggcG93KHB0QS54IC0gcHRCLngsIDIpICsgcG93KHB0QS55IC0gcHRCLnksIDIpICk7XG5cbiAgICAvLyBmaW5kIGxpbmUgZXF1YXRpb24gYXggKyBieSArIGMgPSAwXG4gICAgdmFyIGEgPSAxO1xuICAgIHZhciBiID0gLSAocHRCLnggLSBwdEEueCkgLyAocHRCLnkgLSBwdEEueSk7XG5cbiAgICAvLyBvcmllbnRhdGUgdmVjdG9yIChhLGIpXG4gICAgaWYgKGIgPCAwKXtcbiAgICAgICAgYiA9IC1iO1xuICAgICAgICBhID0gLWE7XG4gICAgfVxuXG4gICAgLy8gbm9ybWFsaXplIHZlY3RvciAoYSxiKVxuICAgIHZhciBuID0gbmV3IFZlY3RvcihhLCBiKTtcbiAgICBuLm5vcm1hbGl6ZSgpO1xuXG4gICAgdmFyIGMgPSAtIChhICogcHRBLnggKyBiICogcHRBLnkpO1xuXG4gICAgLy8gLy8gY2FsY3VsYXRlIHZlY3RvciBkaXJlY3RvclxuICAgIHZhciB2ID0gbmV3IFZlY3RvcihwdEIueCAtIHB0QS54LCBwdEIueSAtIHB0QS55KTtcbiAgICBcbiAgICB2Lm5vcm1hbGl6ZSgpO1xuXG4gICAgdGhpcy5pZCA9IHVuZGVmaW5lZDtcbiAgICB0aGlzLnB0MSA9IHB0QTtcbiAgICB0aGlzLnB0MiA9IHB0QjtcbiAgICB0aGlzLmRpcmVjdGlvbiA9IHY7XG4gICAgdGhpcy5vcnRoRGlyZWN0aW9uID0gbjsgXG4gICAgdGhpcy5kaXN0YW5jZSA9IGRpc3RhbmNlO1xuICAgIHRoaXMucGhlcm9tb24gPSAxL2Rpc3RhbmNlO1xuICAgIHRoaXMubGluZSA9IHtcbiAgICAgICAgYTogYSxcbiAgICAgICAgYjogYixcbiAgICAgICAgYzogYyxcbiAgICB9O1xuXG4gICAgaWYgKHRoaXMuZGlzdGFuY2UgPT09IDApIGNvbnNvbGUubG9nKCdaRVJPICEnKTtcbn1cblxuXG4vLyBzdGF0aWMgbWV0aG9kc1xuRWRnZS5jcmVhdGUgPSBmdW5jdGlvbihwdEEsIHB0Qikge1xuICAgIHZhciBlZGdlID0gbmV3IEVkZ2UocHRBLCBwdEIpO1xuICAgIHJldHVybiBlZGdlO1xufVxuXG5cbi8vIG1ldGhvZHNcbkVkZ2UucHJvdG90eXBlID0ge1xuXG4gICAgZ2V0T3RoZXJQb2ludDogZnVuY3Rpb24ocG9pbnQpIHtcbiAgICAgICAgaWYgKHBvaW50ID09IHRoaXMucHQxKVxuICAgICAgICAgICAgcmV0dXJuIHRoaXMucHQyO1xuICAgICAgICBlbHNlIGlmIChwb2ludCA9PSB0aGlzLnB0MilcbiAgICAgICAgICAgIHJldHVybiB0aGlzLnB0MTtcbiAgICAgICAgZWxzZVxuICAgICAgICAgICAgY29uc29sZS5sb2coXCJFcnJvclwiKTtcbiAgICB9LFxuXG4gICAgY2FsY3VsYXRlRGlzdGFuY2U6IGZ1bmN0aW9uKHgsIHkpIHtcbiAgICAgICAgdmFyIGEgPSB0aGlzLmxpbmUuYSxcbiAgICAgICAgICAgIGIgPSB0aGlzLmxpbmUuYixcbiAgICAgICAgICAgIGMgPSB0aGlzLmxpbmUuYztcbiAgICAgICAgcmV0dXJuIGFicyhhICogeCArIGIgKiB5ICsgYykgLyBNYXRoLnNxcnQoTWF0aC5wb3coYSwyKSArIE1hdGgucG93KGIsMikpO1xuICAgIH0sXG5cbn1cbm1vZHVsZS5leHBvcnRzID0gRWRnZTsiLCIndXNlIHN0cmljdCc7XG5cbnZhciBwYXJzZSA9IHJlcXVpcmUoJ3BhcnNlLXN2Zy1wYXRoJyk7XG5cbnZhciByYW5nZSA9IHJlcXVpcmUoJy4vdXRpbGl0aWVzLmpzJykucmFuZ2U7XG5cbnZhciBQb2ludCA9IHJlcXVpcmUoJy4vcG9pbnQuanMnKTtcbnZhciBzdmdQYXRoID0gcmVxdWlyZSgnLi9zdmdQYXRoLmpzJyk7XG5cbnZhciByYW5kb20gPSBNYXRoLnJhbmRvbTtcblxudmFyIG5iQ2l0eSA9IDI7XG5cbnZhciB0ZXh0TWVzaCA9IHRydWU7XG5cbi8vIEZyYW1lIGRlZmluaXRpb25cbnZhciB4SW5pdCA9IDAsIHlJbml0ID0gMDtcbnZhciB3ID0gMSxcbiAgICBoID0gMTtcblxudmFyIHN2Z1N0cmluZyA9IHN2Z1BhdGg7XG5cbmZ1bmN0aW9uIHN2Z1RvUG9pbnRzKHN2Z1N0cmluZykge1xuICAgIHZhciBwb2ludHMgPSBbXTtcbiAgICB2YXIgZWRnZXMgPSBPYmplY3QuY3JlYXRlKG51bGwpO1xuXG4gICAgdmFyIGJlZ2luaW5nUGF0aDtcblxuICAgIHZhciBYID0gMDtcbiAgICB2YXIgWSA9IDA7XG4gICAgdmFyIG5iUG9pbnRzID0gMDtcbiAgICB2YXIgcHJldlBvaW50O1xuXG4gICAgdmFyIGNvbW1hbmRzID0gcGFyc2Uoc3ZnU3RyaW5nKVxuICAgIGZvciAodmFyIGk9MDsgaTxjb21tYW5kcy5sZW5ndGg7IGkrKyl7XG4gICAgICAgIHZhciBjb21tYW5kID0gY29tbWFuZHNbaV07XG4gICAgICAgIHN3aXRjaCAoY29tbWFuZFswXSkge1xuICAgICAgICAgICAgY2FzZSBcIm1cIjpcbiAgICAgICAgICAgICAgICBYICs9IGNvbW1hbmRbMV07XG4gICAgICAgICAgICAgICAgWSArPSBjb21tYW5kWzJdO1xuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBiZWdpbmluZ1BhdGggPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgIGNhc2UgXCJNXCI6XG4gICAgICAgICAgICAgICAgWCA9IGNvbW1hbmRbMV07XG4gICAgICAgICAgICAgICAgWSA9IGNvbW1hbmRbMl07XG4gICAgICAgICAgICAgICAgcHJldlBvaW50ID0gdW5kZWZpbmVkO1xuICAgICAgICAgICAgICAgIGJlZ2luaW5nUGF0aCA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIGJyZWFrOyBcbiAgICAgICAgICAgIGNhc2UgXCJjXCI6XG4gICAgICAgICAgICAgICAgWCA9IGNvbW1hbmRbMV07XG4gICAgICAgICAgICAgICAgWSA9IGNvbW1hbmRbMl07XG4gICAgICAgICAgICAgICAgcG9pbnRzLnB1c2goe2lkOm5iUG9pbnRzLCB4OlgsIHk6WX0pO1xuXG4gICAgICAgICAgICAgICAgaWYgKHByZXZQb2ludCAhPSB1bmRlZmluZWQpIHtcbiAgICAgICAgICAgICAgICAgICAgZWRnZXNbcHJldlBvaW50XSA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICBuYlBvaW50cysrO1xuICAgICAgICAgICAgICAgIGJyZWFrOyBcbiAgICAgICAgICAgIGNhc2UgXCJsXCI6XG4gICAgICAgICAgICAgICAgWCA9IGNvbW1hbmRbMV07XG4gICAgICAgICAgICAgICAgWSA9IGNvbW1hbmRbMl07XG4gICAgICAgICAgICAgICAgcG9pbnRzLnB1c2goe2lkOm5iUG9pbnRzLCB4OlgsIHk6WX0pO1xuXG4gICAgICAgICAgICAgICAgaWYgKHByZXZQb2ludCAhPSB1bmRlZmluZWQpIHtcbiAgICAgICAgICAgICAgICAgICAgZWRnZXNbcHJldlBvaW50XSA9IG5iUG9pbnRzO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBwcmV2UG9pbnQgPSBuYlBvaW50cztcbiAgICAgICAgICAgICAgICBuYlBvaW50cysrO1xuICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgY2FzZSBcInpcIjpcbiAgICAgICAgICAgICAgICBlZGdlc1twcmV2UG9pbnRdID0gbmJQb2ludHM7XG4gICAgICAgICAgICAgICAgYmVnaW5pbmdQYXRoID0gdW5kZWZpbmVkO1xuICAgICAgICAgICAgICAgIHByZXZQb2ludCA9IHVuZGVmaW5lZDtcbiAgICAgICAgICAgICAgICBicmVhazsgICAgXG4gICAgICAgIH1cbiAgICB9XG4gICAgcmV0dXJuIHtwb2ludHMgOiBwb2ludHMsIGVkZ2VzIDogZWRnZXN9O1xufVxuXG4vLyBpbml0aWFsaXplIHBvaW50c1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKG5iU3RhcnRQb2ludHMsIG5iUmFuZG9tUG9pbnRzKXtcbiAgICB2YXIgcG9pbnRzID0gW107XG4gICAgdmFyIGZvcmNlZEVkZ2VzO1xuICAgIHZhciBjaXR5U2V0O1xuXG4gICAgaWYgKHRleHRNZXNoKXtcblxuICAgICAgICB2YXIgbXlUZXh0ID0gc3ZnVG9Qb2ludHMoc3ZnU3RyaW5nKTtcbiAgICAgICAgcG9pbnRzID0gbXlUZXh0LnBvaW50cztcbiAgICAgICAgZm9yY2VkRWRnZXMgPSBteVRleHQuZWRnZXM7XG4gICAgICAgIGNpdHlTZXQgPSByYW5nZSgwLCBwb2ludHMubGVuZ3RoKTtcblxuICAgICAgICB2YXIgc2NhbGVYID0gMC41O1xuICAgICAgICB2YXIgc2NhbGVZID0gMC40O1xuICAgICAgICB2YXIgZGVsdGFYID0gMC4yNTtcbiAgICAgICAgdmFyIGRlbHRhWSA9IDAuMjU7XG5cbiAgICAgICAgLy8gc2NhbGUgcG9pbnRzIHRvIFswLDFdICsgZGVsdGFcbiAgICAgICAgdmFyIG1heFggPSBNYXRoLm1heC5hcHBseShNYXRoLCBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiBwLnh9KSk7XG4gICAgICAgIHZhciBtaW5YID0gTWF0aC5taW4uYXBwbHkoTWF0aCwgcG9pbnRzLm1hcChmdW5jdGlvbihwKXtyZXR1cm4gcC54fSkpO1xuICAgICAgICB2YXIgbWF4WSA9IE1hdGgubWF4LmFwcGx5KE1hdGgsIHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIHAueX0pKTtcbiAgICAgICAgdmFyIG1pblkgPSBNYXRoLm1pbi5hcHBseShNYXRoLCBwb2ludHMubWFwKGZ1bmN0aW9uKHApe3JldHVybiBwLnl9KSk7XG4gICAgICAgIHBvaW50cyA9IHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7XG4gICAgICAgICAgICB2YXIgeCA9IHNjYWxlWCAqIChwLngtbWluWCkvKG1heFgtbWluWCkgKyBkZWx0YVg7XG4gICAgICAgICAgICB2YXIgeSA9IHNjYWxlWSAqIChwLnktbWluWSkvKG1heFktbWluWSkgKyBkZWx0YVk7XG4gICAgICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgICAgICBuZXdQb2ludC5pZCA9IHAuaWQ7XG5cbiAgICAgICAgICAgIHJldHVybiBuZXdQb2ludDtcbiAgICAgICAgfSk7XG5cbiAgICAgICAgLy8gb25seSBhZGQgcmFuZG9tIHBvaW50c1xuICAgICAgICB2YXIgbmJQb2ludHMgPSBwb2ludHMubGVuZ3RoO1xuICAgICAgICBmb3IodmFyIGk9MDsgaTxuYlJhbmRvbVBvaW50czsgKytpKSB7XG5cbiAgICAgICAgICAgIHZhciB4ID0gcmFuZG9tKCk7XG4gICAgICAgICAgICB2YXIgeSA9IHJhbmRvbSgpO1xuXG4gICAgICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgICAgICBuZXdQb2ludC5pZCA9IG5iUG9pbnRzO1xuXG4gICAgICAgICAgICBwb2ludHMucHVzaChuZXdQb2ludCk7XG5cbiAgICAgICAgICAgIG5iUG9pbnRzKys7XG4gICAgICAgIH1cblxuICAgIH0gZWxzZSB7XG4gICAgICAgIC8vYWRkIHJhbmRvbSBwb2ludHNcblxuICAgICAgICB2YXIgbmJQb2ludHMgPSAwO1xuICAgICAgICBmb3IodmFyIGk9MDsgaTxuYlJhbmRvbVBvaW50czsgKytpKSB7XG5cbiAgICAgICAgICAgIHZhciB4ID0gcmFuZG9tKCk7XG4gICAgICAgICAgICB2YXIgeSA9IHJhbmRvbSgpO1xuXG4gICAgICAgICAgICB2YXIgbmV3UG9pbnQgPSBuZXcgUG9pbnQoeCwgeSk7XG4gICAgICAgICAgICBuZXdQb2ludC5pZCA9IG5iUG9pbnRzO1xuXG4gICAgICAgICAgICBwb2ludHMucHVzaChuZXdQb2ludCk7XG4gICAgICAgICAgICBcbiAgICAgICAgICAgIG5iUG9pbnRzKys7XG4gICAgICAgIH1cblxuICAgICAgICBjaXR5U2V0ID0gcmFuZ2UoMCwgbmJDaXR5KTtcbiAgICAgICAgY29uc29sZS5sb2coY2l0eVNldCk7XG4gICAgfVxuXG5cbiAgICAvLyBpbml0aWFsaXplIHN0YXJ0IHBvaW50c1xuICAgIHZhciBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQgPSBbXTtcblxuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgbmJTdGFydFBvaW50czsgaSsrKXtcbiAgICAgICAgcG9zc2libGVTdGFydFBvaW50c0lkLnB1c2goTWF0aC5mbG9vcihuYlBvaW50cyAqIHJhbmRvbSgpKSk7XG4gICAgfVxuICAgIFxuXG4gICAgcmV0dXJuIHtcbiAgICAgICAgdGV4dE1lc2g6IHRleHRNZXNoLFxuICAgICAgICBwb2ludHM6IHBvaW50cyxcbiAgICAgICAgY2l0eVNldDogY2l0eVNldCxcbiAgICAgICAgcG9zc2libGVTdGFydFBvaW50c0lkOiBwb3NzaWJsZVN0YXJ0UG9pbnRzSWQsXG4gICAgICAgIG5iUmFuZG9tUG9pbnRzOiBuYlJhbmRvbVBvaW50cyxcbiAgICAgICAgZm9yY2VkRWRnZXM6IGZvcmNlZEVkZ2VzXG4gICAgfTtcbn1cbiIsIid1c2Ugc3RyaWN0J1xuXG5tb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uIChjb250YWluZXIpe1xuXG5cdHZhciBtb3VzZSA9IHtcblx0ICAgIHg6IDAsXG5cdCAgICB5OiAwXG5cdH07XG5cblx0Y29udGFpbmVyLmFkZEV2ZW50TGlzdGVuZXIoICdtb3VzZW1vdmUnLCBmdW5jdGlvbihlKXtcblx0ICAgIHZhciByZWN0ID0gY29udGFpbmVyLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpO1xuXG5cdCAgICBtb3VzZS54ID0gKGUuY2xpZW50WCAtIHJlY3QubGVmdCApIC8gcmVjdC53aWR0aDtcblx0ICAgIG1vdXNlLnkgPSAoZS5jbGllbnRZIC0gcmVjdC50b3AgKS8gcmVjdC5oZWlnaHQ7XG5cdH0pO1xuXG5cdHJldHVybiBtb3VzZTtcblxufTtcbiIsIid1c2Ugc3RyaWN0J1xuXG5mdW5jdGlvbiBQb2ludCh4LCB5KSB7XG4gICAgdGhpcy5pZCA9IHVuZGVmaW5lZDsgICAgICAgICAgICAgICAgXG4gICAgdGhpcy54ID0geDtcbiAgICB0aGlzLnkgPSB5O1xuICAgIHRoaXMubmV4dHMgPSBbXTtcbn1cblxubW9kdWxlLmV4cG9ydHMgPSBQb2ludDsiLCIndXNlIHN0cmljdCdcblxudmFyIGFudEZ1bmN0aW9uID0gcmVxdWlyZSgnLi9hbnQuanMnKTtcbnZhciBhbnRzR3JvdXAgPSByZXF1aXJlKCcuL2FudHNHcm91cCcpO1xuXG52YXIgcmFuZG9tID0gTWF0aC5yYW5kb207XG5cbnZhciBSQU5ET01NVlQgPSAwLjAwMztcbnZhciBBTlRTSVpFID0gMC4wMDI7XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oY29udGFpbmVyLCBwb2ludHNNYXAsIG9wdGlvbnMpe1xuXG5cdGlmKCFjb250YWluZXIpXG5cdFx0dGhyb3cgbmV3IFR5cGVFcnJvcignTWlzc2luZyBjb250YWluZXInKTtcblxuXHQvLyBBbnRzIHZhcmlhYmxlc1xuXHR2YXIgZWRnZXMgPSBwb2ludHNNYXAuZWRnZXM7XG5cdHZhciBvYmpQb3B1bGF0aW9uTWF4ID0gb3B0aW9ucy5uYkFudHM7XG5cdHZhciBvYmpQb3B1bGF0aW9uQWRqdXN0ZWQgPSBvYmpQb3B1bGF0aW9uTWF4O1xuXHR2YXIgcG9pbnRzSW5mb3MgPSBwb2ludHNNYXAucG9pbnRzSW5mb3M7XG5cdHZhciBwb3B1bGF0aW9uID0gW107XG5cdHZhciBuYkFudHNQZXJTdGVwID0gMTAwO1xuXHRcblx0dmFyIEFudCA9IGFudEZ1bmN0aW9uKGNvbnRhaW5lciwgcG9pbnRzSW5mb3MsIG9wdGlvbnMpO1xuXHRhbnRzR3JvdXAgPSBhbnRzR3JvdXAoQW50KTtcblxuXHQvLyBBbmltYXRpb24gdmFyaWFibGVzXG5cdHZhciBhbmltSUQ7XG5cdHZhciBkZWx0YVRpbWU7XG5cdHZhciBGUFNDb3VudDtcblx0dmFyIGxhc3RVcGRhdGUgPSBwZXJmb3JtYW5jZS5ub3coKTtcblx0dmFyIEZQU01vbml0b3IgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yKCcjRlBTJyk7XG5cdHZhciBkVE1vbml0b3IgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yKCcjZFQnKTtcblx0dmFyIHdhcm5pbmdNb25pdG9yID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcignI3dhcm5pbmcnKTtcblx0dmFyIHJlZnJlc2hUaW1lID0gMDtcblx0dmFyIG1heERlbHRhVGltZSA9IDQwO1xuXHR2YXIgRlBTT3ZlckxpbWl0Q291bnQgPSAwO1xuXHR2YXIgRlBTVW5kZXJMaW1pdENvdW50ID0gMDtcblxuXHQvLyB3YXJuaW5nIG1lc3NhZ2UgZGlzYXBwZWFycyBhZnRlciA0IHNcblx0d2FybmluZ01vbml0b3IuYWRkRXZlbnRMaXN0ZW5lcihcInRyYW5zaXRpb25lbmRcIiwgZnVuY3Rpb24oKXtcblx0XHR3aW5kb3cuc2V0VGltZW91dChmdW5jdGlvbigpe1xuXHRcdFx0d2FybmluZ01vbml0b3IuY2xhc3NOYW1lID0gXCJpbnZpc2libGVcIjtcblx0XHR9LCA0MDAwKTtcblx0fSlcblxuXHQvLyBDYW52YXNcblx0dmFyIGNhbnZhc0xpc3QgPSBkb2N1bWVudC5nZXRFbGVtZW50c0J5VGFnTmFtZShcImNhbnZhc1wiKTtcblx0XG5cdGlmIChjYW52YXNMaXN0Lmxlbmd0aCA9PT0gMCl7XG5cdFx0dmFyIGNhbnZhcyA9IGRvY3VtZW50LmNyZWF0ZUVsZW1lbnQoXCJjYW52YXNcIik7XG5cdFx0dmFyIHJlY3QgPSBjb250YWluZXIuZ2V0Qm91bmRpbmdDbGllbnRSZWN0KCk7XG5cdFx0Y2FudmFzLndpZHRoID0gcmVjdC53aWR0aDtcblx0XHRjYW52YXMuaGVpZ2h0ID0gcmVjdC5oZWlnaHQ7XG5cdFx0Y2FudmFzLnN0eWxlLmJhY2tncm91bmRDb2xvciA9IFwicmdiYSgyNTAsIDI1MCwgMjUwLCAwKVwiOyBcblx0XHRjb250YWluZXIuYXBwZW5kQ2hpbGQoY2FudmFzKTtcblx0fVxuXHRlbHNle1xuXHRcdHZhciBjYW52YXMgPSBjYW52YXNMaXN0WzBdO1xuXHRcdGNvbnNvbGUubG9nKCdDQU5WQVMnKTtcblx0fVxuXHR2YXIgY29udGV4dCA9IGNhbnZhcy5nZXRDb250ZXh0KFwiMmRcIik7XG5cdGNvbnRleHQuY2xlYXJSZWN0ICggMCAsIDAgLCBjYW52YXMud2lkdGgsIGNhbnZhcy5oZWlnaHQgKTtcblx0XG5cblx0ZnVuY3Rpb24gY2hlY2tBbnROdW1iZXIoYW50TnVtYmVyKXtcblx0XHRpZiAoYW50TnVtYmVyIDwgb2JqUG9wdWxhdGlvbkFkanVzdGVkIC0gNTApe1xuXHRcdFx0RlBTTW9uaXRvci5zdHlsZS5jb2xvciA9IFwiZ3JlZW5cIjtcblx0XHRcdHBvcHVsYXRpb24gPSBhbnRzR3JvdXAuY3JlYXRlKHBvcHVsYXRpb24pO1xuXHRcdH1cdFxuXHRcdGVsc2UgaWYgKGFudE51bWJlciA+IG9ialBvcHVsYXRpb25BZGp1c3RlZCl7XG5cdFx0XHRwb3B1bGF0aW9uID0gYW50c0dyb3VwLnJlbW92ZShwb3B1bGF0aW9uLCBhbnROdW1iZXIgLSBvYmpQb3B1bGF0aW9uQWRqdXN0ZWQpO1xuXHRcdFx0RlBTTW9uaXRvci5zdHlsZS5jb2xvciA9IFwicmVkXCI7XG5cdFx0XHR3YXJuaW5nTW9uaXRvci5jbGFzc05hbWUgPSBcInZpc2libGVcIjtcblx0XHR9XG5cdFx0ZWxzZXtcblx0XHRcdEZQU01vbml0b3Iuc3R5bGUuY29sb3IgPSBcIndoaXRlXCI7XG5cdFx0XHQvLyB3YXJuaW5nTW9uaXRvci5jbGFzc05hbWUgPSBcImludmlzaWJsZVwiO1xuXHRcdH1cblx0fVxuXG5cdGZ1bmN0aW9uIGRpc3BsYXlGUFMoZFQpe1xuXHRcdEZQU0NvdW50ID0gKDEwMDAvZFQpLnRvRml4ZWQoMik7XG5cdFx0dmFyIHQgPSBkVC50b0ZpeGVkKDIpO1xuXHRcdEZQU01vbml0b3IudGV4dENvbnRlbnQgPSAnRlBTIDogJyArIEZQU0NvdW50OyAgXG5cdFx0ZFRNb25pdG9yLnRleHRDb250ZW50ID0gJ25iQW50cyA6ICcgKyBwb3B1bGF0aW9uLmxlbmd0aDtcblx0XHQvLyBkVE1vbml0b3IuaW5uZXJUZXh0ID0gJ2RUIDogJyArIHQgKyAnbXMnO1xuXHR9XG5cblx0ZnVuY3Rpb24gdGljaygpIHtcblx0XHR2YXIgbm93ID0gcGVyZm9ybWFuY2Uubm93KCk7XG5cdFx0ZGVsdGFUaW1lID0gbm93IC0gbGFzdFVwZGF0ZTtcblx0XHRsYXN0VXBkYXRlID0gbm93O1xuXHRcdHJlZnJlc2hUaW1lICs9IGRlbHRhVGltZS8xMDAwOyAvLyBpbiBzZWNvbmRzXG5cblx0XHQvLyBjb25zb2xlLmxvZygnbmJBbnRzJywgcG9wdWxhdGlvbi5sZW5ndGgpO1xuXG5cdFx0Y2hlY2tBbnROdW1iZXIocG9wdWxhdGlvbi5sZW5ndGgpO1xuXG5cdFx0Ly8gZGlzcGxheSBGUFMgaW5mbyBldmVyeSAwLjMgc1xuXHRcdGlmIChyZWZyZXNoVGltZSA+IDAuMyl7XG5cdFx0XHRkaXNwbGF5RlBTKGRlbHRhVGltZSk7XG5cdFx0XHRyZWZyZXNoVGltZSA9IDA7IFxuXHRcdH1cblxuXHRcdC8vIHJlbW92ZSBhbnRzIHdoZW4gZnJhbWUgcmF0ZSBpcyB0b28gbG93XG5cdFx0aWYgKEZQU092ZXJMaW1pdENvdW50ID09PSAxMCkge1xuXHRcdFx0b2JqUG9wdWxhdGlvbkFkanVzdGVkID0gb2JqUG9wdWxhdGlvbkFkanVzdGVkICogbWF4RGVsdGFUaW1lIC8gZGVsdGFUaW1lO1xuXHRcdFx0RlBTT3ZlckxpbWl0Q291bnQgPSAwO1xuXHRcdH1cblxuXHRcdHdoaWxlIChGUFNVbmRlckxpbWl0Q291bnQgPiA1MCAmJiBvYmpQb3B1bGF0aW9uQWRqdXN0ZWQgPCBvYmpQb3B1bGF0aW9uTWF4KSB7XG5cdFx0XHRvYmpQb3B1bGF0aW9uQWRqdXN0ZWQgKz0gMTA7XG5cdFx0fVxuXG5cdFx0Ly8gY2hlY2sgZHVyYXRpb24gb2Ygb3Zlci91bmRlciBmcmFtZXJhdGUgbGltaXQgcGVyaW9kc1xuXHRcdGlmIChkZWx0YVRpbWUgPiBtYXhEZWx0YVRpbWUpe1xuXHRcdFx0RlBTT3ZlckxpbWl0Q291bnQrKztcblx0XHRcdEZQU1VuZGVyTGltaXRDb3VudCA9IDA7XG5cdFx0fVxuXHRcdGVsc2Uge1xuXHRcdFx0RlBTT3ZlckxpbWl0Q291bnQgPSAwO1xuXHRcdFx0RlBTVW5kZXJMaW1pdENvdW50Kys7XG5cdFx0fVxuXG5cdFx0Ly8gZHJhdyBpbiBjYW52YXNcblx0XHR2YXIgdyA9IGNhbnZhcy53aWR0aDtcblx0XHR2YXIgaCA9IGNhbnZhcy5oZWlnaHQ7XG5cdFx0dmFyIG1vdXNlID0gW2xhc3RNb3VzZU1vdmVFdmVudC5jbGllbnRYL3csIGxhc3RNb3VzZU1vdmVFdmVudC5jbGllbnRZL2hdO1xuXHRcdGNvbnRleHQuc2V0VHJhbnNmb3JtKHcsIDAsIDAsIGgsIDAsIDApO1xuXHRcdGNvbnRleHQuZmlsbFN0eWxlID0gXCJyZ2JhKDI1MCwgMjUwLCAyNTAsIDAuNClcIjtcblx0XHRjb250ZXh0LmZpbGxSZWN0KDAsMCx3LGgpO1xuXG5cdFx0Ly8gLy8gZWRnZXNcblx0XHQvLyBjb250ZXh0LnN0cm9rZVN0eWxlID0gXCIjMDAwXCI7XG5cdFx0Ly8gZm9yKHZhciBpPTA7IGkgPCBlZGdlcy5sZW5ndGg7ICsraSkge1xuXHRcdC8vICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IDAuMDAwMTtcblx0XHQvLyAgICAgdmFyIGVkZ2UgPSBlZGdlc1tpXTtcblx0XHQvLyAgICAgaWYgKGVkZ2UucGhlcm9tb24gIT0gMCl7XG5cdFx0Ly8gICAgICAgICBjb250ZXh0LmxpbmVXaWR0aCA9IE1hdGgubWluKDAuMDAwMDEgKiBlZGdlLnBoZXJvbW9uLCAwLjAxKTtcblx0XHQvLyAgICAgfSBlbHNlIHtcblx0XHQvLyAgICAgICAgIGNvbnRleHQubGluZVdpZHRoID0gMC4wMDAwMTtcblx0XHQvLyAgICAgfVxuXHRcdC8vICAgICBjb250ZXh0LmJlZ2luUGF0aCgpO1xuXHRcdC8vICAgICBjb250ZXh0Lm1vdmVUbyhwb2ludHNJbmZvcy5wb2ludHNbZWRnZS5wdDEuaWRdLngsIHBvaW50c0luZm9zLnBvaW50c1tlZGdlLnB0MS5pZF0ueSk7XG5cdFx0Ly8gICAgIGNvbnRleHQubGluZVRvKHBvaW50c0luZm9zLnBvaW50c1tlZGdlLnB0Mi5pZF0ueCwgcG9pbnRzSW5mb3MucG9pbnRzW2VkZ2UucHQyLmlkXS55KTtcblx0XHQvLyAgICAgY29udGV4dC5zdHJva2UoKTtcblx0XHQvLyB9XG5cblx0XHQvLyAvLyB2ZXJ0aWNlc1xuXHRcdC8vIGZvcih2YXIgaT0wOyBpPHBvaW50c0luZm9zLnBvaW50cy5sZW5ndGg7ICsraSkge1xuXHRcdC8vICAgICBjb250ZXh0LmJlZ2luUGF0aCgpXG5cdFx0Ly8gICAgIHZhciBwb2ludCA9IHBvaW50c0luZm9zLnBvaW50c1tpXTtcblx0XHQvLyAgICAgaWYgKHBvaW50c0luZm9zLmNpdHlTZXQuaW5kZXhPZihwb2ludC5pZCkgIT0gLTEpe1xuXHRcdC8vICAgICAgICAgY29udGV4dC5maWxsU3R5bGUgPSBcIiMwMTAxREZcIjtcblx0XHQvLyAgICAgICAgIGNvbnRleHQuYXJjKHBvaW50LngsIHBvaW50LnksIDAuMDA2LCAwLCAyKk1hdGguUEkpO1xuXHRcdC8vICAgICB9XG5cdFx0Ly8gICAgIGVsc2Uge1xuXHRcdC8vICAgICAgICAgY29udGV4dC5maWxsU3R5bGUgPSBcIiMwMDBcIjtcblx0XHQvLyAgICAgICAgIGNvbnRleHQuYXJjKHBvaW50c0luZm9zLnBvaW50c1tpXS54LCBwb2ludHNJbmZvcy5wb2ludHNbaV0ueSwgMC4wMDMsIDAsIDIqTWF0aC5QSSk7XG5cdFx0Ly8gICAgIH1cblx0XHQvLyAgICAgY29udGV4dC5jbG9zZVBhdGgoKTtcblx0XHQvLyAgICAgY29udGV4dC5maWxsKCk7XG5cdFx0Ly8gfVxuXG5cdFx0Ly8gbW92ZSBhbnRzXG5cdFx0cG9wdWxhdGlvbi5mb3JFYWNoKGZ1bmN0aW9uKGFudCl7XG5cdFx0XHRhbnQudHJhbnNpdCgpO1xuXHRcdH0pO1xuXG5cdFx0Ly8gcGhlcm9tb24gZXZhcG9yYXRpb25cblx0XHRlZGdlcy5mb3JFYWNoKGZ1bmN0aW9uKGVkZ2Upe1xuXHRcdFx0aWYoZWRnZS5waGVyb21vbiA+IDApe1xuXHRcdFx0XHRlZGdlLnBoZXJvbW9uIC09IDAuMDAwMTtcblx0XHRcdH1cblx0XHR9KTtcblxuXHRcdC8vIGFudHNcblx0XHRwb3B1bGF0aW9uLmZvckVhY2goZnVuY3Rpb24oYW50KXtcblx0XHRcdGNvbnRleHQuYmVnaW5QYXRoKClcblx0XHRcdHZhciB4ID0gYW50LnggKyBSQU5ET01NVlQqcmFuZG9tKCk7XG5cdFx0XHR2YXIgeSA9IGFudC55ICsgUkFORE9NTVZUKnJhbmRvbSgpO1xuXG5cdFx0XHRjb250ZXh0LmZpbGxTdHlsZSA9IFwiYmxhY2tcIlxuXHRcdFx0Y29udGV4dC5maWxsUmVjdCh4LCB5LCBBTlRTSVpFLCBBTlRTSVpFKTtcblx0XHRcdGNvbnRleHQuY2xvc2VQYXRoKCk7XG5cdFx0XHRjb250ZXh0LmZpbGwoKTtcblx0XHR9KVxuXHR9O1xuXHRcblx0dmFyIGxhc3RNb3VzZU1vdmVFdmVudCA9IHtcblx0XHRjbGllbnRYOiAwLFxuXHRcdGNsaWVudFk6IDBcblx0fTtcblx0XG5cdGNvbnRhaW5lci5hZGRFdmVudExpc3RlbmVyKCdtb3VzZW1vdmUnLCBmdW5jdGlvbihlKXtcblx0XHRsYXN0TW91c2VNb3ZlRXZlbnQgPSBlO1xuXHR9KTtcblx0XG5cdHZhciBwYXVzZWQgPSBmYWxzZTtcblx0XG5cdGZ1bmN0aW9uIHRvZ2dsZVBsYXlQYXVzZSgpe1xuXHRcdHBhdXNlZCA9ICFwYXVzZWQ7XG5cdFx0aWYoIXBhdXNlZClcblx0XHRcdGFuaW1hdGUoKTtcblx0fVxuXG5cdGZ1bmN0aW9uIHJlc2V0KCl7XG5cdFx0cG9wdWxhdGlvbiA9IFtdO1xuXHRcdGVkZ2VzID0gW107XG5cdFx0cG9pbnRzSW5mb3MgPSBbXTtcblxuXHRcdGNhbmNlbEFuaW1hdGlvbkZyYW1lKGFuaW1JRCk7XG5cdH1cblx0XG5cdC8vIGNvbnRhaW5lci5hZGRFdmVudExpc3RlbmVyKCdjbGljaycsIHRvZ2dsZVBsYXlQYXVzZSk7XG5cblx0ZnVuY3Rpb24gYW5pbWF0ZSgpe1xuXHRcdHRpY2soKTtcblx0XHRcblx0XHRpZighcGF1c2VkKVxuXHRcdFx0YW5pbUlEID0gcmVxdWVzdEFuaW1hdGlvbkZyYW1lKGFuaW1hdGUpO1xuXHR9XG5cdGFuaW1hdGUoKTtcblxuXHRmdW5jdGlvbiBtb2RpZnlBbnRzKG9wdHMpe1xuXHRcdG9ialBvcHVsYXRpb25NYXggPSBvcHRzLm5iQW50cztcblx0XHRvYmpQb3B1bGF0aW9uQWRqdXN0ZWQgPSBvYmpQb3B1bGF0aW9uTWF4O1xuXG5cdFx0cG9wdWxhdGlvbi5mb3JFYWNoKGZ1bmN0aW9uKGFudCl7XG5cdFx0XHRhbnQudmVsb2NpdHkgPSBvcHRzLnZlbG9jaXR5O1xuXHRcdFx0YW50LmludGVsbGlnZW5jZSA9IG9wdHMuaW50ZWxsaWdlbmNlO1xuXHRcdFx0YW50LnJlcFNpemUgPSBvcHRzLnJlcFNpemU7XG5cdFx0XHRhbnQucmVwU3BlZWQgPSBvcHRzLnJlcFNwZWVkO1xuXHRcdH0pO1xuXHR9XG5cdFxuXHRyZXR1cm4ge1xuXHRcdHRvZ2dsZVBsYXlQYXVzZTogdG9nZ2xlUGxheVBhdXNlLFxuXHRcdHJlc2V0OiByZXNldCxcblx0XHQvLyBzaG91bGQgYmUgYSBnZXR0ZXIvc2V0dGVyLCBidXQgSUU4XG5cdFx0Z2V0QW50Q291bnQ6IGZ1bmN0aW9uKCl7XG5cdFx0XHRyZXR1cm4gcG9wdWxhdGlvbi5sZW5ndGg7XG5cdFx0fSxcblx0XHRtb2RpZnlBbnRzOiBtb2RpZnlBbnRzXG5cdH1cbn1cbiIsIid1c2Ugc3RyaWN0JztcblxucmVxdWlyZSgnLi4vaW5kZXguanMnKShkb2N1bWVudC5ib2R5LnF1ZXJ5U2VsZWN0b3IoJy5jb2xvbnknKSk7IiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgcG9pbnQxID0gXCJtIDE4LjI1LDE5LjUgYyAxOC4yNSwxNS43Njc3Njg4IDE1LjAwMDA5NTEsMTIuNzUgMTEsMTIuNzUgYyA2Ljk5OTkwNDg4LDEyLjc1IDMuNzUsMTUuNzY3NzY4OCAzLjc1LDE5LjUgYyAzLjc1LDIzLjIzMjIzMTIgNi45OTk5MDQ4OCwyNi4yNSAxMSwyNi4yNSBjIDE1LjAwMDA5NTEsMjYuMjUgMTguMjUsMjMuMjMyMjMxMiAxOC4yNSwxOS41IHogbSA0LjI1LDE5LjUgYyA0LjI1LDE2LjA1MjUyOTQgNy4yNjgxMDg2MywxMy4yNSAxMSwxMy4yNSBjIDE0LjczMTg5MTQsMTMuMjUgMTcuNzUsMTYuMDUyNTI5NCAxNy43NSwxOS41IGMgMTcuNzUsMjIuOTQ3NDcwNiAxNC43MzE4OTE0LDI1Ljc1IDExLDI1Ljc1IGMgNy4yNjgxMDg2MywyNS43NSA0LjI1LDIyLjk0NzQ3MDYgNC4yNSwxOS41IHpcIjtcbnZhciBwb2ludDIgPSBcIm0gODkuMjUsOC41IGMgODkuMjUsNC4yMTYwNTc4NyA4NS41NTI4NzE2LDAuNzUgODEsMC43NSBjIDc2LjQ0NzEyODQsMC43NSA3Mi43NSw0LjIxNjA1Nzg3IDcyLjc1LDguNSBjIDcyLjc1LDEyLjc4Mzk0MjEgNzYuNDQ3MTI4NCwxNi4yNSA4MSwxNi4yNSBjIDg1LjU1Mjg3MTYsMTYuMjUgODkuMjUsMTIuNzgzOTQyMSA4OS4yNSw4LjUgeiBtIDczLjI1LDguNSBjIDczLjI1LDQuNDk5NjcwODggNzYuNzE2MzE1NiwxLjI1IDgxLDEuMjUgYyA4NS4yODM2ODQ0LDEuMjUgODguNzUsNC40OTk2NzA4OCA4OC43NSw4LjUgYyA4OC43NSwxMi41MDAzMjkxIDg1LjI4MzY4NDQsMTUuNzUgODEsMTUuNzUgYyA3Ni43MTYzMTU2LDE1Ljc1IDczLjI1LDEyLjUwMDMyOTEgNzMuMjUsOC41IHpcIjtcbnZhciBwb2ludDMgPSBcIm0gMTYwLjI1LDExIGMgMTYwLjI1LDUuMzMzNDgzNzUgMTU1LjIwODE2OCwwLjc1IDE0OSwwLjc1IGMgMTQyLjc5MTgzMiwwLjc1IDEzNy43NSw1LjMzMzQ4Mzc1IDEzNy43NSwxMSBjIDEzNy43NSwxNi42NjY1MTYyIDE0Mi43OTE4MzIsMjEuMjUgMTQ5LDIxLjI1IGMgMTU1LjIwODE2OCwyMS4yNSAxNjAuMjUsMTYuNjY2NTE2MiAxNjAuMjUsMTEgeiBtIDEzOC4yNSwxMSBjIDEzOC4yNSw1LjYyMDgyMTI1IDE0My4wNTc5MDMsMS4yNSAxNDksMS4yNSBjIDE1NC45NDIwOTcsMS4yNSAxNTkuNzUsNS42MjA4MjEyNSAxNTkuNzUsMTEgYyAxNTkuNzUsMTYuMzc5MTc4NyAxNTQuOTQyMDk3LDIwLjc1IDE0OSwyMC43NSBjIDE0My4wNTc5MDMsMjAuNzUgMTM4LjI1LDE2LjM3OTE3ODcgMTM4LjI1LDExIHpcIjtcbnZhciBwb2ludDQgPSBcIm0gMTYwLjI1LDc2LjUgYyAxNjAuMjUsNzIuNzY3NzY4OCAxNTcuMDAwMDk1LDY5Ljc1IDE1Myw2OS43NSBjIDE0OC45OTk5MDUsNjkuNzUgMTQ1Ljc1LDcyLjc2Nzc2ODggMTQ1Ljc1LDc2LjUgYyAxNDUuNzUsODAuMjMyMjMxMiAxNDguOTk5OTA1LDgzLjI1IDE1Myw4My4yNSBjIDE1Ny4wMDAwOTUsODMuMjUgMTYwLjI1LDgwLjIzMjIzMTIgMTYwLjI1LDc2LjUgeiBtIDE0Ni4yNSw3Ni41IGMgMTQ2LjI1LDczLjA1MjUyOTQgMTQ5LjI2ODEwOSw3MC4yNSAxNTMsNzAuMjUgYyAxNTYuNzMxODkxLDcwLjI1IDE1OS43NSw3My4wNTI1Mjk0IDE1OS43NSw3Ni41IGMgMTU5Ljc1LDc5Ljk0NzQ3MDYgMTU2LjczMTg5MSw4Mi43NSAxNTMsODIuNzUgYyAxNDkuMjY4MTA5LDgyLjc1IDE0Ni4yNSw3OS45NDc0NzA2IDE0Ni4yNSw3Ni41IHpcIjtcbnZhciBwb2ludDUgPSBcIm0gOTUuMjUsNzYgYyA5NS4yNSw3MC4zMzYyNzk1IDkwLjQzNDQwNjUsNjUuNzUgODQuNSw2NS43NSBjIDc4LjU2NTU5MzUsNjUuNzUgNzMuNzUsNzAuMzM2Mjc5NSA3My43NSw3NiBjIDczLjc1LDgxLjY2MzcyMDUgNzguNTY1NTkzNSw4Ni4yNSA4NC41LDg2LjI1IGMgOTAuNDM0NDA2NSw4Ni4yNSA5NS4yNSw4MS42NjM3MjA1IDk1LjI1LDc2IHogbSA3NC4yNSw3NiBjIDc0LjI1LDcwLjYxODAyNTUgNzguODM2NDI2OCw2Ni4yNSA4NC41LDY2LjI1IGMgOTAuMTYzNTczMiw2Ni4yNSA5NC43NSw3MC42MTgwMjU1IDk0Ljc1LDc2IGMgOTQuNzUsODEuMzgxOTc0NSA5MC4xNjM1NzMyLDg1Ljc1IDg0LjUsODUuNzUgYyA3OC44MzY0MjY4LDg1Ljc1IDc0LjI1LDgxLjM4MTk3NDUgNzQuMjUsNzYgelwiO1xudmFyIHBvaW50NiA9IFwibSAyMC4yNSw3NSBjIDIwLjI1LDcwLjk5MTkzMzggMTYuNzc2NDk5NSw2Ny43NSAxMi41LDY3Ljc1IGMgOC4yMjM1MDA0Niw2Ny43NSA0Ljc1LDcwLjk5MTkzMzggNC43NSw3NSBjIDQuNzUsNzkuMDA4MDY2MiA4LjIyMzUwMDQ2LDgyLjI1IDEyLjUsODIuMjUgYyAxNi43NzY0OTk1LDgyLjI1IDIwLjI1LDc5LjAwODA2NjIgMjAuMjUsNzUgeiBtIDUuMjUsNzUgYyA1LjI1LDcxLjI3NjA3OTcgOC40OTIyMjgyOSw2OC4yNSAxMi41LDY4LjI1IGMgMTYuNTA3NzcxNyw2OC4yNSAxOS43NSw3MS4yNzYwNzk3IDE5Ljc1LDc1IGMgMTkuNzUsNzguNzIzOTIwMyAxNi41MDc3NzE3LDgxLjc1IDEyLjUsODEuNzUgYyA4LjQ5MjIyODI5LDgxLjc1IDUuMjUsNzguNzIzOTIwMyA1LjI1LDc1IHpcIjtcbnZhciBwb2ludDcgPSBcIm0gMjMuMjUsMTM5IGMgMjMuMjUsMTMyLjc4Njc5NyAxOC4yMTMyMDM0LDEyNy43NSAxMiwxMjcuNzUgYyA1Ljc4Njc5NjU2LDEyNy43NSAwLjc1LDEzMi43ODY3OTcgMC43NSwxMzkgYyAwLjc1LDE0NS4yMTMyMDMgNS43ODY3OTY1NiwxNTAuMjUgMTIsMTUwLjI1IGMgMTguMjEzMjAzNCwxNTAuMjUgMjMuMjUsMTQ1LjIxMzIwMyAyMy4yNSwxMzkgeiBtIDEuMjUsMTM5IGMgMS4yNSwxMzMuMDYyOTM5IDYuMDYyOTM4OTQsMTI4LjI1IDEyLDEyOC4yNSBjIDE3LjkzNzA2MTEsMTI4LjI1IDIyLjc1LDEzMy4wNjI5MzkgMjIuNzUsMTM5IGMgMjIuNzUsMTQ0LjkzNzA2MSAxNy45MzcwNjExLDE0OS43NSAxMiwxNDkuNzUgYyA2LjA2MjkzODk0LDE0OS43NSAxLjI1LDE0NC45MzcwNjEgMS4yNSwxMzkgelwiO1xudmFyIHBvaW50OCA9IFwibSA5NS4yNDI5MjA5LDEzMy4xNDE2MjQgYyA5NS41NTM3OTI3LDEyNy4yMDk4MzYgOTAuNzU2MjY5OCwxMjIuMTQxNzM2IDg0LjUzMzc0MDQsMTIxLjgxNTYyNyBjIDc4LjMxMTIxMSwxMjEuNDg5NTE4IDczLjAxMDIwODYsMTI2LjAyODM3NyA3Mi42OTkzMzY4LDEzMS45NjAxNjUgYyA3Mi4zODg0NjUsMTM3Ljg5MTk1MyA3Ny4xODU5ODc5LDE0Mi45NjAwNTMgODMuNDA4NTE3MywxNDMuMjg2MTYyIGMgODkuNjMxMDQ2NywxNDMuNjEyMjcxIDk0LjkzMjA0OSwxMzkuMDczNDEyIDk1LjI0MjkyMDksMTMzLjE0MTYyNCB6IG0gNzMuMTk4NjUxNiwxMzEuOTg2MzMzIGMgNzMuNDk0NzcxMSwxMjYuMzM2MDM2IDc4LjU1NTM4OCwxMjIuMDAzMDAxIDg0LjUwNzU3MjQsMTIyLjMxNDk0MiBjIDkwLjQ1OTc1NjcsMTIyLjYyNjg4MyA5NS4wMzk3MjU2LDEyNy40NjUxNTkgOTQuNzQzNjA2MSwxMzMuMTE1NDU2IGMgOTQuNDQ3NDg2NiwxMzguNzY1NzU0IDg5LjM4Njg2OTYsMTQzLjA5ODc4OCA4My40MzQ2ODUzLDE0Mi43ODY4NDcgYyA3Ny40ODI1MDA5LDE0Mi40NzQ5MDcgNzIuOTAyNTMyLDEzNy42MzY2MyA3My4xOTg2NTE2LDEzMS45ODYzMzMgelwiO1xudmFyIHBvaW50OSA9IFwibSAxNjcuNzI4MDIzLDEzNS4wNjI0NzkgYyAxNjguMDM4MzI3LDEyOS4xNDE1MzMgMTYzLjQ4MjExNywxMjQuMDg5ODY0IDE1Ny41NTE2NjIsMTIzLjc3OTA2MiBjIDE1MS42MjEyMDgsMTIzLjQ2ODI2MSAxNDYuNTYxOTE0LDEyOC4wMTYwMDIgMTQ2LjI1MTYxLDEzMy45MzY5NDggYyAxNDUuOTQxMzA3LDEzOS44NTc4OTQgMTUwLjQ5NzUxNywxNDQuOTA5NTYyIDE1Ni40Mjc5NzEsMTQ1LjIyMDM2NCBjIDE2Mi4zNTg0MjYsMTQ1LjUzMTE2NiAxNjcuNDE3NzIsMTQwLjk4MzQyNSAxNjcuNzI4MDIzLDEzNS4wNjI0NzkgeiBtIDE0Ni43NTA5MjUsMTMzLjk2MzExNiBjIDE0Ny4wNDY3NjcsMTI4LjMxODEyMSAxNTEuODcwNjE3LDEyMy45ODIwMTggMTU3LjUyNTQ5NCwxMjQuMjc4Mzc3IGMgMTYzLjE4MDM3MiwxMjQuNTc0NzM3IDE2Ny41MjQ1NSwxMjkuMzkxMzE3IDE2Ny4yMjg3MDksMTM1LjAzNjMxMSBjIDE2Ni45MzI4NjcsMTQwLjY4MTMwNSAxNjIuMTA5MDE3LDE0NS4wMTc0MDkgMTU2LjQ1NDEzOSwxNDQuNzIxMDUgYyAxNTAuNzk5MjYyLDE0NC40MjQ2OSAxNDYuNDU1MDgzLDEzOS42MDgxMSAxNDYuNzUwOTI1LDEzMy45NjMxMTYgelwiO1xudmFyIGxldHRyZXNBTlQgPSBcIm0gMjUuNjE4NTQxNCwzMS45NDQ3MjY2IGwgMjUuNjE4NTQxNCwzMS4xMDAzOTA2IGwgMjUuMzg4MjY4LDMwLjQwOTU3MDMgbCAyNC45Mjc3MjExLDMwLjAyNTc4MTIgbCAyNC41NDM5MzIxLDI5LjU2NTIzNDQgbCAyNC4wODMzODUyLDI5LjE4MTQ0NTMgbCAyMy42OTk1OTYxLDI4LjcyMDg5ODQgbCAyMy4wMDg3NzU4LDI4LjQ5MDYyNSBsIDIyLjE0NTM2NzgsMjguNDkwNjI1IGwgMjEuMTUxMjcxMSwyOC40OTA2MjUgbCAyMC40NzU3NjgsMjguNDkwNjI1IGwgMTkuODYxNzA1NSwyOC43MjA4OTg0IGwgMTkuNDc3OTE2NCwyOS4xODE0NDUzIGwgMTguNzg3MDk2MSwyOS4zMzQ5NjA5IGwgMTcuNTQ3ODU3OSwyOS4zMzQ5NjA5IGwgMTYuNjY5OTA3LDI5LjMzNDk2MDkgbCAxNi4wMDI0OTcxLDI5LjMzNDk2MDkgbCAxNS4zMzQ5ODAyLDI5LjMzNDk2MDkgbCAxNC41NjU0MTY0LDI5LjMzNDk2MDkgbCAxMy45NTEzNTM5LDI5LjU2NTIzNDQgbCAxMy41Njc1NjQ5LDMwLjAyNTc4MTIgbCAxMi44NzY3NDQ2LDMwLjE3OTI5NjkgbCAxMi4wMzI0MDg2LDMwLjE3OTI5NjkgbCAxMS40MTgzNDYxLDMwLjQwOTU3MDMgbCAxMS4wMzQ1NTcxLDMwLjg3MDExNzIgbCAxMC4zNDM3MzY3LDMxLjAyMzYzMjggbCA5LjcyOTY3NDI0LDMxLjI1MzkwNjIgbCA5LjM0NTg4NTE4LDMxLjcxNDQ1MzEgbCA4Ljg4NTMzODMsMzIuMDk4MjQyMiBsIDguNTAxNTQ5MjQsMzIuNTU4Nzg5MSBsIDguMDQxMDAyMzcsMzIuOTQyNTc4MSBsIDcuODg3NDg2NzQsMzMuNjMzMzk4NCBsIDcuNjU3MjEzMywzNC4yNDc0NjA5IGwgNy4xOTY2NjY0MywzNC42MzEyNSBsIDcuMDQzMTUwOCwzNS4zMjIwNzAzIGwgNy4wNDMxNTA4LDM2LjE2NjQwNjIgbCA3LjE5NjY2NjQzLDM2Ljc4MDQ2ODcgbCA3LjY1NzIxMzMsMzcuMTY0MjU3OCBsIDguMDQxMDAyMzcsMzcuNjI0ODA0NyBsIDguNTAxNTQ5MjQsMzguMDA4NTkzNyBsIDguODg1MzM4MywzOC40NjkxNDA2IGwgOS4zNDU4ODUxOCwzOC44NTI5Mjk3IGwgOS43Mjk2NzQyNCwzOS4zMTM0NzY2IGwgMTAuMzQzNzM2NywzOS40NjY5OTIyIGwgMTEuNDcyNTk2MSwzOS40NjY5OTIyIGwgMTIuMzEyMjk2MSwzOS40NjY5OTIyIGwgMTIuODc0MzgxOCwzOS40NjY5OTIyIGwgMTMuNjQ5Njg3MywzOS40NjY5OTIyIGwgMTQuNTM5ODUyOCwzOS40NjY5OTIyIGwgMTUuNDA5NzUyNCwzOS40NjY5OTIyIGwgMTYuMTAwNTcyNywzOS4zMTM0NzY2IGwgMTYuNDg0MzYxNywzOC44NTI5Mjk3IGwgMTcuMDk4NDI0MiwzOC42MjI2NTYyIGwgMTcuNzg5MjQ0NiwzOC40NjkxNDA2IGwgMTguMTczMDMzNiwzOC4wMDg1OTM3IGwgMTguNjMzNTgwNSwzNy42MjQ4MDQ3IGwgMTkuMDE3MzY5NiwzNy4xNjQyNTc4IGwgMTkuNjMxNDMyMSwzNi45MzM5ODQ0IGwgMjAuMzIyMjUyNCwzNi43ODA0Njg3IGwgMjAuNzA2MDQxNCwzNi4zMTk5MjE5IGwgMjEuMzIwMTAzOSwzNi4wODk2NDg0IGwgMjIuMDEwOTI0MiwzNS45MzYxMzI4IGwgMjIuMzk0NzEzMywzNS40NzU1ODU5IGwgMjMuMDA4Nzc1OCwzNS4yNDUzMTI1IGwgMjMuNjk5NTk2MSwzNS4wOTE3OTY5IGwgMjMuOTI5ODY5NiwzNC40Nzc3MzQ0IGwgMjQuMDgzMzg1MiwzMy43ODY5MTQxIGwgMjQuNTQzOTMyMSwzMy40MDMxMjUgbCAyNC45Mjc3MjExLDMyLjk0MjU3ODEgbCAyNS4zODgyNjgsMzIuNTU4Nzg5MSBsIDI1LjYxODU0MTQsMzEuOTQ0NzI2NiB6IG0gMzcuNTkyNzYwMiw0My45MTg5NDUzIGwgMzcuMjA4OTcxMSw0NC4zNzk0OTIyIGwgMzYuNzQ4NDI0Miw0NC43NjMyODEzIGwgMzYuMzY0NjM1Miw0NS4yMjM4MjgxIGwgMzUuNjczODE0OSw0NS4zNzczNDM4IGwgMzUuMDU5NzUyNCw0NS42MDc2MTcyIGwgMzQuNjc1OTYzMyw0Ni4wNjgxNjQxIGwgMzMuOTg1MTQzLDQ2LjIyMTY3OTcgbCAzMy4zOTUxNzIsNDYuMjIxNjc5NyBsIDMyLjQ4NzQyMTksNDYuMjIxNjc5NyBsIDMxLjQ1MjEzNTIsNDYuMjIxNjc5NyBsIDMwLjgzODA3MjcsNDYuMDY4MTY0MSBsIDMwLjQ1NDI4MzYsNDUuNjA3NjE3MiBsIDI5Ljc2MzQ2MzMsNDUuMzc3MzQzOCBsIDI5LjE0OTQwMDgsNDUuMjIzODI4MSBsIDI4Ljc2NTYxMTcsNDQuNzYzMjgxMyBsIDI4LjMwNTA2NDksNDQuMzc5NDkyMiBsIDI3LjkyMTI3NTgsNDMuOTE4OTQ1MyBsIDI3LjIzMDQ1NTUsNDMuNjg4NjcxOSBsIDI2LjYxNjM5Myw0My41MzUxNTYzIGwgMjYuMjMyNjAzOSw0My4wNzQ2MDk0IGwgMjUuNzcyMDU3MSw0My4wNzQ2MDk0IGwgMjUuMzg4MjY4LDQzLjUzNTE1NjMgbCAyNC42OTc0NDc3LDQzLjY4ODY3MTkgbCAyMy44NTMxMTE3LDQzLjY4ODY3MTkgbCAyMy4yMzkwNDkyLDQzLjkxODk0NTMgbCAyMi44NTUyNjAyLDQ0LjM3OTQ5MjIgbCAyMi4xNjQ0Mzk5LDQ0LjUzMzAwNzggbCAyMS41NTAzNzc0LDQ0Ljc2MzI4MTMgbCAyMS4xNjY1ODgzLDQ1LjIyMzgyODEgbCAyMC40NzU3NjgsNDUuMzc3MzQzOCBsIDE5LjYzMTQzMjEsNDUuMzc3MzQzOCBsIDE5LjAxNzM2OTYsNDUuNjA3NjE3MiBsIDE4LjYzMzU4MDUsNDYuMDY4MTY0MSBsIDE3Ljk0Mjc2MDIsNDYuMjIxNjc5NyBsIDE3LjEyNzkwODcsNDYuMjIxNjc5NyBsIDE2LjI1NDA4ODMsNDYuMjIxNjc5NyBsIDE1LjY0MDAyNTgsNDYuNDUxOTUzMSBsIDE1LjI1NjIzNjcsNDYuOTEyNSBsIDE0Ljc5NTY4OTksNDYuOTEyNSBsIDE0LjQxMTkwMDgsNDYuNDUxOTUzMSBsIDEzLjcyMTA4MDUsNDYuMjIxNjc5NyBsIDEyLjc5MjY4MzcsNDYuMjIxNjc5NyBsIDEyLjAxNjk4LDQ2LjIyMTY3OTcgbCAxMS4zNjQ2OSw0Ni4yMjE2Nzk3IGwgMTAuNDg0OTEsNDYuMjIxNjc5NyBsIDkuNDk5NDAwOCw0Ni4yMjE2Nzk3IGwgOC44ODUzMzgzLDQ2LjA2ODE2NDEgbCA4LjUwMTU0OTI0LDQ1LjYwNzYxNzIgbCA3LjgxMDcyODkzLDQ1LjM3NzM0MzggbCA3LjE5NjY2NjQzLDQ1LjIyMzgyODEgbCA2LjgxMjg3NzM3LDQ0Ljc2MzI4MTMgbCA2LjEyMjA1NzA1LDQ0LjUzMzAwNzggbCA1LjUwNzk5NDU1LDQ0LjM3OTQ5MjIgbCA1LjEyNDIwNTQ5LDQzLjkxODk0NTMgbCA0LjY2MzY1ODYyLDQzLjUzNTE1NjMgbCA0LjI3OTg2OTU1LDQzLjA3NDYwOTQgbCAzLjgxOTMyMjY4LDQyLjY5MDgyMDMgbCAzLjQzNTUzMzYyLDQyLjIzMDI3MzQgbCAyLjk3NDk4Njc0LDQxLjg0NjQ4NDQgbCAyLjU5MTE5NzY4LDQxLjM4NTkzNzUgbCAyLjEzMDY1MDgsNDEuMDAyMTQ4NCBsIDEuOTc3MTM1MTgsNDAuMzg4MDg1OSBsIDEuNzQ2ODYxNzQsMzkuNjk3MjY1NiBsIDEuMjg2MzE0ODcsMzkuMzEzNDc2NiBsIDEuMTMyNzk5MjQsMzguNjk5NDE0MSBsIDEuMTMyNzk5MjQsMzcuODU1MDc4MSBsIDAuOTAyNTI1ODA0LDM3LjE2NDI1NzggbCAwLjQ0MTk3ODkyOSwzNi43ODA0Njg4IGwgMC4yODg0NjMzMDQsMzYuMTY2NDA2MyBsIDAuMjg4NDYzMzA0LDM1LjI1MjE0MzMgbCAwLjI4ODQ2MzMwNCwzNC4zOTcwODc1IGwgMC4yODg0NjMzMDQsMzMuNjMzMzk4NCBjIDAuMjg4NDYzMzA0LDMzLjUwMDg3MDUgMC40NDE5Nzg5MjksMzIuOTQyNTc4MSAwLjQ0MTk3ODkyOSwzMi45NDI1NzgxIGwgMC45MDI1MjU4MDQsMzIuNTU4Nzg5MSBsIDEuMTMyNzk5MjQsMzEuOTQ0NzI2NiBsIDEuMTMyNzk5MjQsMzEuMTAwMzkwNiBsIDEuMjg2MzE0ODcsMzAuNDA5NTcwMyBsIDEuNzQ2ODYxNzQsMzAuMDI1NzgxMyBsIDEuOTc3MTM1MTgsMjkuNDExNzE4OCBsIDIuMTMwNjUwOCwyOC43MjA4OTg0IGwgMi41OTExOTc2OCwyOC4zMzcxMDk0IGwgMi45NzQ5ODY3NCwyNy44NzY1NjI1IGwgMy40MzU1MzM2MiwyNy40OTI3NzM0IGwgMy44MTkzMjI2OCwyNy4wMzIyMjY2IGwgNC4yNzk4Njk1NSwyNi42NDg0Mzc1IGwgNC42NjM2NTg2MiwyNi4xODc4OTA2IGwgNS4yNzc3MjExMiwyNS45NTc2MTcyIGwgNi4xMjIwNTcwNSwyNS45NTc2MTcyIGwgNi44MTI4NzczNywyNS44MDQxMDE2IGwgNy4xOTY2NjY0MywyNS4zNDM1NTQ3IGwgNy44MTA3Mjg5MywyNS4xMTMyODEzIGwgOC42NTUwNjQ4NywyNS4xMTMyODEzIGwgOS4zNDU4ODUxOCwyNC45NTk3NjU2IGwgOS43Mjk2NzQyNCwyNC40OTkyMTg4IGwgMTAuMzQzNzM2NywyNC4yNjg5NDUzIGwgMTAuOTM3MDU4NSwyNC4yNjg5NDUzIGwgMTEuNjY5NTI5NCwyNC4yNjg5NDUzIGwgMTIuMTI0NDQwNCwyNC4yNjg5NDUzIGwgMTIuNzg5MDkyNSwyNC4yNjg5NDUzIGwgMTMuNzIxMDgwNSwyNC4yNjg5NDUzIGMgMTMuOTU2NjA4NCwyNC4yNjg5NDUzIDE0LjQxMTkwMDgsMjQuMTE1NDI5NyAxNC40MTE5MDA4LDI0LjExNTQyOTcgbCAxNC43OTU2ODk5LDIzLjY1NDg4MjggbCAxNS40MDk3NTI0LDIzLjQyNDYwOTQgbCAxNS44OTE0OTI5LDIzLjQyNDYwOTQgbCAxNi4zNzgwMzQyLDIzLjQyNDYwOTQgbCAxNi43MTY1MDY4LDIzLjQyNDYwOTQgbCAxNy4zMjY0NjE4LDIzLjQyNDYwOTQgbCAxNy45NDI3NjAyLDIzLjQyNDYwOTQgYyAxOC4xNzgyODgxLDIzLjQyNDYwOTQgMTguNjMzNTgwNSwyMy4yNzEwOTM4IDE4LjYzMzU4MDUsMjMuMjcxMDkzOCBsIDE5LjAxNzM2OTYsMjIuODEwNTQ2OSBsIDE5LjYzMTQzMjEsMjIuNTgwMjczNCBsIDIwLjEyNTg4NzUsMjIuNTgwMjczNCBsIDIwLjY1NTY2MTYsMjIuNTgwMjczNCBsIDIxLjMwMjYzMzgsMjIuNTgwMjczNCBsIDIyLjQ5MTY3MzIsMjIuNTgwMjczNCBsIDIzLjQzNDc4MzEsMjIuNTgwMjczNCBsIDIzLjg1MzExMTcsMjIuNTgwMjczNCBsIDI0LjU0MzkzMjEsMjIuNDI2NzU3OCBsIDI0LjkyNzcyMTEsMjEuOTY2MjEwOSBsIDI1LjM4ODI2OCwyMS41ODI0MjE5IGwgMjUuNzcyMDU3MSwyMS4xMjE4NzUgbCAyNi4zODYxMTk2LDIwLjg5MTYwMTYgbCAyNy4wNzY5Mzk5LDIwLjczODA4NTkgbCAyNy4zMDcyMTMzLDIwLjEyNDAyMzQgbCAyNy4zMDcyMTMzLDE5LjMyMDEwNjcgbCAyNy4zMDcyMTMzLDE4LjQzNTM1MTYgbCAyNy4wNzY5Mzk5LDE3Ljc0NDUzMTMgbCAyNi42MTYzOTMsMTcuMzYwNzQyMiBsIDI2LjYxNjM5MywxNi45MDAxOTUzIGwgMjcuMDc2OTM5OSwxNi41MTY0MDYzIGwgMjcuMDc2OTM5OSwxNi4wNTU4NTk0IGwgMjYuNjE2MzkzLDE1LjY3MjA3MDMgbCAyNi40NjI4Nzc0LDE1LjA1ODAwNzggbCAyNi4yMzI2MDM5LDE0LjM2NzE4NzUgbCAyNS43NzIwNTcxLDEzLjk4MzM5ODQgbCAyNS4zODgyNjgsMTMuNTIyODUxNiBsIDI0LjkyNzcyMTEsMTMuMTM5MDYyNSBsIDI0LjU0MzkzMjEsMTIuNjc4NTE1NiBsIDIzLjg1MzExMTcsMTIuNDQ4MjQyMiBsIDIzLjAwODc3NTgsMTIuNDQ4MjQyMiBsIDIyLjM5NDcxMzMsMTIuMjk0NzI2NiBsIDIyLjAxMDkyNDIsMTEuODM0MTc5NyBsIDIxLjU1MDM3NzQsMTEuODM0MTc5NyBsIDIxLjE2NjU4ODMsMTIuMjk0NzI2NiBsIDIwLjQ3NTc2OCwxMi40NDgyNDIyIGwgMTkuNTM1NjE1MSwxMi40NDgyNDIyIGwgMTguNjc4OTAxOCwxMi40NDgyNDIyIGwgMTguMTAwNjQ2MSwxMi40NDgyNDIyIGwgMTcuNDkwNjkxMSwxMi40NDgyNDIyIGwgMTYuODY5MjAyNywxMi40NDgyNDIyIGwgMTYuMjk2MjY0OSwxMi40NDgyNDIyIGwgMTUuNzU1ODU1MiwxMi40NDgyNDIyIGwgMTUuMTY0NDA4OCwxMi40NDgyNDIyIGwgMTQuNTY1NDE2NCwxMi40NDgyNDIyIGwgMTMuOTUxMzUzOSwxMi42Nzg1MTU2IGwgMTMuNTY3NTY0OSwxMy4xMzkwNjI1IGwgMTMuMTA3MDE4LDEzLjUyMjg1MTYgbCAxMi43MjMyMjg5LDEzLjk4MzM5ODQgbCAxMi4yNjI2ODIxLDE0LjM2NzE4NzUgbCAxMi4xMDkxNjY0LDE1LjA1ODAwNzggbCAxMi4xMDkxNjY0LDE1LjcwMjA4OTkgbCAxMi4xMDkxNjY0LDE2LjY5NDU3OSBsIDEyLjEwOTE2NjQsMTcuNTkxMDE1NiBsIDExLjg3ODg5MywxOC4yMDUwNzgxIGwgMTEuNDE4MzQ2MSwxOC41ODg4NjcyIGwgMTEuMDM0NTU3MSwxOS4wNDk0MTQxIGwgMTAuNTc0MDEwMiwxOS40MzMyMDMxIGwgMTAuMTkwMjIxMSwxOS44OTM3NSBsIDkuNDk5NDAwOCwyMC4wNDcyNjU2IGwgOC42NDQ0Nzg0LDIwLjA0NzI2NTYgbCA3LjgxMDcyODkzLDIwLjA0NzI2NTYgbCA3LjE5NjY2NjQzLDE5Ljg5Mzc1IGwgNi44MTI4NzczNywxOS40MzMyMDMxIGwgNi4xMjIwNTcwNSwxOS4yMDI5Mjk3IGwgNS41MDc5OTQ1NSwxOS4wNDk0MTQxIGwgNS4xMjQyMDU0OSwxOC41ODg4NjcyIGwgNC42NjM2NTg2MiwxOC4yMDUwNzgxIGwgNC41MTAxNDI5OSwxNy41OTEwMTU2IGwgNC41MTAxNDI5OSwxNi43NDY2Nzk3IGwgNC4yNzk4Njk1NSwxNi4wNTU4NTk0IGwgMy44MTkzMjI2OCwxNS42NzIwNzAzIGwgMy42NjU4MDcwNSwxNS4wNTgwMDc4IGwgMy42NjU4MDcwNSwxNC4yMTM2NzE5IGwgMy44MTkzMjI2OCwxMy41MjI4NTE2IGwgNC4yNzk4Njk1NSwxMy4xMzkwNjI1IGwgNC41MTAxNDI5OSwxMi41MjUgbCA0LjY2MzY1ODYyLDExLjgzNDE3OTcgbCA1LjEyNDIwNTQ5LDExLjQ1MDM5MDYgbCA1LjUwNzk5NDU1LDEwLjk4OTg0MzggbCA1Ljk2ODU0MTQzLDEwLjYwNjA1NDcgbCA2LjM1MjMzMDQ5LDEwLjE0NTUwNzggbCA2Ljk2NjM5Mjk5LDkuOTE1MjM0MzggbCA3LjY1NzIxMzMsOS43NjE3MTg3NSBsIDguMDQxMDAyMzcsOS4zMDExNzE4OCBsIDguNjU1MDY0ODcsOS4wNzA4OTg0NCBsIDkuMzQ1ODg1MTgsOC45MTczODI4MSBsIDkuNzI5Njc0MjQsOC40NTY4MzU5NCBsIDEwLjM0MzczNjcsOC4yMjY1NjI1IGwgMTAuOTg5MDYyMSw4LjIyNjU2MjUgbCAxMS40MzYxMDAxLDguMjI2NTYyNSBsIDExLjk0MjgwNzUsOC4yMjY1NjI1IGwgMTIuMjY4MTU4Myw4LjIyNjU2MjUgbCAxMi45NTQ5NzkzLDguMjI2NTYyNSBsIDEzLjcyMTA4MDUsOC4yMjY1NjI1IGwgMTQuNDExOTAwOCw4LjA3MzA0Njg4IGwgMTQuNzk1Njg5OSw3LjYxMjUgbCAxNS4yNTYyMzY3LDcuNjEyNSBsIDE1LjY0MDAyNTgsOC4wNzMwNDY4OCBsIDE2LjI1NDA4ODMsOC4yMjY1NjI1IGwgMTYuOTYwOTE3LDguMjI2NTYyNSBsIDE3LjQ5MTUxOTksOC4yMjY1NjI1IGwgMTguMDI5MDk4LDguMjI2NTYyNSBsIDE4Ljc5MzI2ODMsOC4yMjY1NjI1IGwgMTkuNjMxNDMyMSw4LjIyNjU2MjUgbCAyMC4zMjIyNTI0LDguNDU2ODM1OTQgbCAyMC43MDYwNDE0LDguOTE3MzgyODEgbCAyMS4zMjAxMDM5LDkuMDcwODk4NDQgbCAyMi4wNDQwODI4LDkuMDcwODk4NDQgbCAyMi44Nzk0NTYsOS4wNzA4OTg0NCBsIDIzLjg1MzExMTcsOS4wNzA4OTg0NCBjIDI0LjA5NTA0NzQsOS4wNzA4OTg0NCAyNC41NDM5MzIxLDkuMzAxMTcxODggMjQuNTQzOTMyMSw5LjMwMTE3MTg4IGwgMjQuOTI3NzIxMSw5Ljc2MTcxODc1IGwgMjUuNTQxNzgzNiw5LjkxNTIzNDM4IGwgMjYuMjMyNjAzOSwxMC4xNDU1MDc4IGwgMjYuNjE2MzkzLDEwLjYwNjA1NDcgbCAyNy4wNzY5Mzk5LDEwLjk4OTg0MzggbCAyNy40NjA3Mjg5LDExLjQ1MDM5MDYgbCAyNy45MjEyNzU4LDExLjgzNDE3OTcgbCAyOC4zMDUwNjQ5LDEyLjI5NDcyNjYgbCAyOC43NjU2MTE3LDEyLjY3ODUxNTYgbCAyOS4xNDk0MDA4LDEzLjEzOTA2MjUgbCAyOS42MDk5NDc3LDEzLjUyMjg1MTYgbCAyOS44NDAyMjExLDE0LjIxMzY3MTkgbCAyOS45OTM3MzY3LDE0LjgyNzczNDQgbCAzMC40NTQyODM2LDE1LjIxMTUyMzQgbCAzMC42ODQ1NTcxLDE1LjkwMjM0MzggbCAzMC42ODQ1NTcxLDE2LjUwMzI3NzQgYyAzMC42ODQ1NTcxLDE2LjcyNzE1MzIgMzAuNjg0NTU3MSwxNi45NTEwMjkgMzAuNjg0NTU3MSwxNy4xNzQ5MDQ4IGwgMzAuNjg0NTU3MSwxNy44NzY1NzQxIGwgMzAuNjg0NTU3MSwxOC41NDIwNTQ5IGwgMzAuNjg0NTU3MSwxOS4xNjUyMDA4IGMgMzAuNjg0NTU3MSwxOS43MDkxMTQ3IDMwLjY4NDU1NzEsMjAuMTI0MDIzNCAzMC42ODQ1NTcxLDIwLjEyNDAyMzQgbCAzMC44MzgwNzI3LDIwLjczODA4NTkgbCAzMS4yOTg2MTk2LDIxLjEyMTg3NSBsIDMxLjUyODg5MywyMS44MTI2OTUzIGwgMzEuNTI4ODkzLDIyLjYwNjA3MzcgYyAzMS41Mjg4OTMsMjIuODYyODY5IDMxLjUyODg5MywyMy4xMTk2NjQzIDMxLjUyODg5MywyMy4zNzY0NTk2IGwgMzEuNTI4ODkzLDI0LjExMDU4OCBjIDMxLjUyODg5MywyNC4zNzcwNzQyIDMxLjUyODg5MywyNC42NTU1ODQ2IDMxLjUyODg5MywyNC45NDM0NzUgbCAzMS41Mjg4OTMsMjUuODg3NDEzNiBsIDMxLjUyODg5MywyNi43NzQxNjg5IGMgMzEuNTI4ODkzLDI3LjA3NDQwNDEgMzEuNTI4ODkzLDI3LjM3NDYzOTQgMzEuNTI4ODkzLDI3LjY3NDg3NDcgbCAzMS41Mjg4OTMsMjguNTg3MTEzOCBsIDMxLjUyODg5MywyOS4yNjIxMjUyIGwgMzEuNTI4ODkzLDMwLjAxOTI1MTIgYyAzMS41Mjg4OTMsMzAuMjY0NDc2OCAzMS41Mjg4OTMsMzAuNTEyODA2MyAzMS41Mjg4OTMsMzAuNzYxMjUyNyBsIDMxLjUyODg5MywzMS41OTE2NTM0IGwgMzEuNTI4ODkzLDMyLjI3NTY0MjkgbCAzMS41Mjg4OTMsMzMuMTUwMDM2MSBsIDMxLjUyODg5MywzMy44MTM4NTk0IGwgMzEuNTI4ODkzLDM0LjU1MzM3NDYgYyAzMS41Mjg4OTMsMzUuMDMxMDc3OCAzMS41Mjg4OTMsMzUuMzIyMDcwMyAzMS41Mjg4OTMsMzUuMzIyMDcwMyBsIDMxLjY4MjQwODYsMzUuOTM2MTMyOCBsIDMyLjE0Mjk1NTUsMzYuMzE5OTIxOSBsIDMyLjUyNjc0NDYsMzYuNzgwNDY4OCBsIDMyLjk4NzI5MTQsMzcuMTY0MjU3OCBsIDMzLjM3MTA4MDUsMzcuNjI0ODA0NyBsIDMzLjk4NTE0MywzNy43NzgzMjAzIGwgMzQuNjc1OTYzMywzOC4wMDg1OTM4IGwgMzUuMDU5NzUyNCwzOC40NjkxNDA2IGwgMzUuNjczODE0OSwzOC42MjI2NTYzIGwgMzYuMzY0NjM1MiwzOC44NTI5Mjk3IGwgMzYuNzQ4NDI0MiwzOS4zMTM0NzY2IGwgMzcuMzYyNDg2NywzOS40NjY5OTIyIGwgMzguMDUzMzA3MSwzOS42OTcyNjU2IGwgMzguMjgzNTgwNSw0MC4zODgwODU5IGwgMzguMjgzNTgwNSw0MS4wMjEyMzIzIGwgMzguMjgzNTgwNSw0Mi4wNzUzOTM3IGwgMzguMjgzNTgwNSw0Mi45MjEwOTM4IGwgMzguMDUzMzA3MSw0My41MzUxNTYzIGwgMzcuNTkyNzYwMiw0My45MTg5NDUzIHogbSA4Ny40MTgxODU3LDQwLjU0MTYwMTYgbCA4Ny4yNjQ2NzAxLDQxLjIzMjQyMTkgbCA4Ny4wMzQzOTY2LDQxLjg0NjQ4NDQgbCA4Ni41NzM4NDk3LDQyLjIzMDI3MzQgbCA4Ni4xOTAwNjA3LDQyLjY5MDgyMDMgbCA4NS43Mjk1MTM4LDQzLjA3NDYwOTQgbCA4NS4zNDU3MjQ3LDQzLjUzNTE1NjMgbCA4NC42NTQ5MDQ0LDQzLjY4ODY3MTkgbCA4NC4wNDA4NDE5LDQzLjkxODk0NTMgbCA4My42NTcwNTI5LDQ0LjM3OTQ5MjIgbCA4Mi45NjYyMzI2LDQ0LjUzMzAwNzggbCA4Mi4zODI4MTkxLDQ0LjUzMzAwNzggbCA4MS43NjE1NDYzLDQ0LjUzMzAwNzggbCA4MS4xNDkzOTkzLDQ0LjUzMzAwNzggbCA4MC4zOTkxMTIyLDQ0LjUzMzAwNzggbCA3OS43Njk4NzY3LDQ0LjUzMzAwNzggbCA3OS4yNjc0MTg4LDQ0LjUzMzAwNzggbCA3OC41OTkyNjQyLDQ0LjUzMzAwNzggbCA3Ny45NDcxMjQ2LDQ0LjUzMzAwNzggbCA3Ny4yNTk3MzQyLDQ0LjUzMzAwNzggbCA3Ni42Njk0MTc3LDQ0LjUzMzAwNzggbCA3Ni4xNTg5OTcsNDQuNTMzMDA3OCBsIDc1LjU4NDYwNiw0NC41MzMwMDc4IGwgNzUuMDM0MTkyNiw0NC41MzMwMDc4IGwgNzQuNTA3NzU3LDQ0LjUzMzAwNzggbCA3NC4yNzEyMDA5LDQ0LjUzMzAwNzggbCA3My43MDk2MDM5LDQ0LjUzMzAwNzggbCA3Mi44MzQyMDEzLDQ0LjUzMzAwNzggYyA3Mi44Njg2NDEsNDQuNTAwMDQxMyA3Mi4yMjAxMzg4LDQ0LjM3OTQ5MjIgNzIuMjIwMTM4OCw0NC4zNzk0OTIyIGwgNzEuODM2MzQ5Nyw0My45MTg5NDUzIGwgNzEuMzc1ODAyOSw0My41MzUxNTYzIGwgNzEuMjIyMjg3Miw0Mi45MjEwOTM4IGwgNzAuOTkyMDEzOCw0Mi4yMzAyNzM0IGwgNzAuNTMxNDY2OSw0MS44NDY0ODQ0IGwgNzAuMzc3OTUxMyw0MS4yMzI0MjE5IGwgNzAuNTMxNDY2OSw0MC41NDE2MDE2IGwgNzAuOTkyMDEzOCw0MC4xNTc4MTI1IGwgNzEuMzc1ODAyOSwzOS42OTcyNjU2IGwgNzEuODM2MzQ5NywzOS4zMTM0NzY2IGwgNzIuMjIwMTM4OCwzOC44NTI5Mjk3IGwgNzIuNjgwNjg1NywzOC40NjkxNDA2IGwgNzMuMDY0NDc0NywzOC4wMDg1OTM4IGwgNzMuNjc4NTM3MiwzNy43NzgzMjAzIGwgNzQuMzY5MzU3NiwzNy42MjQ4MDQ3IGwgNzQuNzUzMTQ2NiwzNy4xNjQyNTc4IGwgNzUuMjEzNjkzNSwzNi43ODA0Njg4IGwgNzUuNTk3NDgyNiwzNi4zMTk5MjE5IGwgNzYuMDU4MDI5NCwzNS45MzYxMzI4IGwgNzYuMjg4MzAyOSwzNS4zMjIwNzAzIGwgNzYuMjg4MzAyOSwzNC40Nzc3MzQ0IGwgNzYuNDQxODE4NSwzMy43ODY5MTQxIGwgNzYuOTAyMzY1NCwzMy40MDMxMjUgbCA3Ny4xMzI2Mzg4LDMyLjc4OTA2MjUgbCA3Ny4xMzI2Mzg4LDMyLjA1ODY1MTcgbCA3Ny4xMzI2Mzg4LDMxLjEzNDc5NDcgbCA3Ny4xMzI2Mzg4LDMwLjA2ODA1NTcgbCA3Ny4xMzI2Mzg4LDI5LjI2MzAxMzYgbCA3Ny4xMzI2Mzg4LDI4LjM4NzAyMjUgbCA3Ny4xMzI2Mzg4LDI3LjUwNDk0NzUgbCA3Ny4xMzI2Mzg4LDI2LjYwMjgzMTUgbCA3Ny4xMzI2Mzg4LDI1Ljc3MTY2NDQgbCA3Ny4xMzI2Mzg4LDI0LjgyMzc0MDMgbCA3Ny4xMzI2Mzg4LDIzLjkyMTcxMzggbCA3Ny4xMzI2Mzg4LDIzLjA5MDQ1NzIgbCA3Ny4xMzI2Mzg4LDIyLjEzMTM0OTQgbCA3Ny4xMzI2Mzg4LDIxLjE2NDI3ODkgbCA3Ny4xMzI2Mzg4LDIwLjIzODI3NDcgbCA3Ny4xMzI2Mzg4LDE5LjU0MjkyMTUgbCA3Ny4xMzI2Mzg4LDE4LjkyMTY0ODggbCA3Ny4xMzI2Mzg4LDE4LjE3OTQxMzkgbCA3Ny4xMzI2Mzg4LDE3LjU5MTAxNTYgbCA3Ny4yODYxNTQ0LDE2LjkwMDE5NTMgbCA3Ny43NDY3MDEzLDE2LjUxNjQwNjMgbCA3Ny45NzY5NzQ3LDE1LjkwMjM0MzggbCA3Ny43NDY3MDEzLDE1LjIxMTUyMzQgbCA3Ny4yODYxNTQ0LDE0LjgyNzczNDQgbCA3Ny4xMzI2Mzg4LDE0LjIxMzY3MTkgbCA3Ni45MDIzNjU0LDEzLjUyMjg1MTYgbCA3Ni40NDE4MTg1LDEzLjEzOTA2MjUgbCA3Ni4wNTgwMjk0LDEyLjY3ODUxNTYgbCA3NS4zNjcyMDkxLDEyLjQ0ODI0MjIgbCA3NC43NTMxNDY2LDEyLjI5NDcyNjYgbCA3NC4zNjkzNTc2LDExLjgzNDE3OTcgbCA3My42Nzg1MzcyLDExLjYwMzkwNjMgbCA3My4wNjQ0NzQ3LDExLjQ1MDM5MDYgbCA3Mi42ODA2ODU3LDEwLjk4OTg0MzggbCA3MS45ODk4NjU0LDEwLjc1OTU3MDMgbCA3MS4yMTM1MDgxLDEwLjc1OTU3MDMgbCA3MC4yODg0ODgsMTAuNzU5NTcwMyBsIDY5Ljc3MTE3ODIsMTAuNzU5NTcwMyBsIDY5LjA1Mzk5NDYsMTAuNzU5NTcwMyBsIDY4LjI5NjgxODQsMTAuNzU5NTcwMyBsIDY3LjU2NTg1NjYsMTAuNzU5NTcwMyBsIDY2LjkwNDU5MTIsMTAuNzU5NTcwMyBsIDY2LjA3OTUxMzgsMTAuNzU5NTcwMyBjIDY2LjExNjM4OTMsMTAuODEwMjAxIDY1LjQ2NTQ1MTMsMTAuOTg5ODQzOCA2NS40NjU0NTEzLDEwLjk4OTg0MzggbCA2NS4wODE2NjIyLDExLjQ1MDM5MDYgbCA2NC4zOTA4NDE5LDExLjYwMzkwNjMgbCA2My43NzY3Nzk0LDExLjgzNDE3OTcgbCA2My42MjMyNjM4LDEyLjUyNSBsIDYzLjM5Mjk5MDQsMTMuMTM5MDYyNSBsIDYyLjkzMjQ0MzUsMTMuNTIyODUxNiBsIDYyLjU0ODY1NDQsMTMuOTgzMzk4NCBsIDYyLjA4ODEwNzYsMTQuMzY3MTg3NSBsIDYxLjcwNDMxODUsMTQuODI3NzM0NCBsIDYxLjI0Mzc3MTYsMTUuMjExNTIzNCBsIDYxLjA5MDI1NiwxNS45MDIzNDM4IGwgNjAuODU5OTgyNiwxNi41MTY0MDYzIGwgNjAuMzk5NDM1NywxNi45MDAxOTUzIGwgNjAuMDE1NjQ2NiwxNy4zNjA3NDIyIGwgNTkuNTU1MDk5NywxNy43NDQ1MzEzIGwgNTkuNDAxNTg0MSwxOC40MzUzNTE2IGwgNTkuMTcxMzEwNywxOS4wNDk0MTQxIGwgNTguNzEwNzYzOCwxOS40MzMyMDMxIGwgNTguNTU3MjQ4MiwyMC4xMjQwMjM0IGwgNTguMzI2OTc0NywyMC43MzgwODU5IGwgNTcuODY2NDI3OSwyMS4xMjE4NzUgbCA1Ny43MTI5MTIyLDIxLjgxMjY5NTMgbCA1Ny40ODI2Mzg4LDIyLjQyNjc1NzggbCA1Ny4wMjIwOTE5LDIyLjgxMDU0NjkgbCA1Ni44Njg1NzYzLDIzLjUwMTM2NzIgbCA1Ni44Njg1NzYzLDIzLjk4OTcxMDIgbCA1Ni44Njg1NzYzLDI0LjY5ODkzMSBsIDU2Ljg2ODU3NjMsMjUuNjU0ODE3OSBsIDU2Ljg2ODU3NjMsMjYuMzI1MTE5NyBsIDU2Ljg2ODU3NjMsMjcuMDkwMjU4NyBsIDU2Ljg2ODU3NjMsMjcuODAwNDYzNyBsIDU2Ljg2ODU3NjMsMjguNTY3MzgyOCBsIDU3LjAyMjA5MTksMjkuMTgxNDQ1MyBsIDU3LjQ4MjYzODgsMjkuNTY1MjM0NCBsIDU3LjcxMjkxMjIsMzAuMjU2MDU0NyBsIDU3LjcxMjkxMjIsMzEuMTAwMzkwNiBsIDU3LjQ4MjYzODgsMzEuNzE0NDUzMSBsIDU3LjAyMjA5MTksMzIuMDk4MjQyMiBsIDU2Ljg2ODU3NjMsMzIuNzg5MDYyNSBsIDU3LjAyMjA5MTksMzMuNDAzMTI1IGwgNTcuNDgyNjM4OCwzMy43ODY5MTQxIGwgNTcuNzEyOTEyMiwzNC40Nzc3MzQ0IGwgNTcuNzEyOTEyMiwzNS4zMjIwNzAzIGwgNTcuODY2NDI3OSwzNS45MzYxMzI4IGwgNTguMzI2OTc0NywzNi4zMTk5MjE5IGwgNTguNzEwNzYzOCwzNi43ODA0Njg4IGwgNTkuMTcxMzEwNywzNy4xNjQyNTc4IGwgNTkuNTU1MDk5NywzNy42MjQ4MDQ3IGwgNjAuMTY5MTYyMiwzNy43NzgzMjAzIGwgNjAuODU5OTgyNiwzOC4wMDg1OTM4IGwgNjEuMjQzNzcxNiwzOC40NjkxNDA2IGwgNjEuNzA0MzE4NSwzOC44NTI5Mjk3IGwgNjIuMDg4MTA3NiwzOS4zMTM0NzY2IGwgNjIuNTQ4NjU0NCwzOS42OTcyNjU2IGwgNjIuNzc4OTI3OSw0MC4zODgwODU5IGwgNjIuOTMyNDQzNSw0MS4wMDIxNDg0IGwgNjMuMzkyOTkwNCw0MS4zODU5Mzc1IGwgNjMuMzkyOTkwNCw0MS44NDY0ODQ0IGwgNjIuOTMyNDQzNSw0Mi4yMzAyNzM0IGwgNjIuNzc4OTI3OSw0Mi45MjEwOTM4IGwgNjIuNTQ4NjU0NCw0My41MzUxNTYzIGwgNjEuODU3ODM0MSw0My42ODg2NzE5IGwgNjEuMjQzNzcxNiw0My45MTg5NDUzIGwgNjAuODU5OTgyNiw0NC4zNzk0OTIyIGwgNjAuMTY5MTYyMiw0NC41MzMwMDc4IGwgNTkuMzI0ODI2Myw0NC41MzMwMDc4IGwgNTguNzEwNzYzOCw0NC43NjMyODEzIGwgNTguMzI2OTc0Nyw0NS4yMjM4MjgxIGwgNTcuODY2NDI3OSw0NS4yMjM4MjgxIGwgNTcuNDgyNjM4OCw0NC43NjMyODEzIGwgNTcuMDIyMDkxOSw0NC43NjMyODEzIGwgNTYuNjM4MzAyOSw0NS4yMjM4MjgxIGwgNTUuOTQ3NDgyNiw0NS4zNzczNDM4IGwgNTUuMTAzMTQ2Niw0NS4zNzczNDM4IGwgNTQuNDg5MDg0MSw0NS4yMjM4MjgxIGwgNTQuMTA1Mjk1MSw0NC43NjMyODEzIGwgNTMuNjQ0NzQ4Miw0NC43NjMyODEzIGwgNTMuMjYwOTU5MSw0NS4yMjM4MjgxIGwgNTIuNTcwMTM4OCw0NS4zNzczNDM4IGwgNTEuNDk3NDkxMSw0NS4zNzczNDM4IGwgNTAuNTkwNjMzMiw0NS4zNzczNDM4IGwgNDkuOTg1Mzc1NCw0NS4zNzczNDM4IGwgNDkuMTkyNzk1MSw0NS4zNzczNDM4IGMgNDkuMjI2ODQwMiw0NS4zNDQzMDM5IDQ4LjU3ODczMjYsNDUuMjIzODI4MSA0OC41Nzg3MzI2LDQ1LjIyMzgyODEgbCA0OC4xOTQ5NDM1LDQ0Ljc2MzI4MTMgbCA0Ny41MDQxMjMyLDQ0LjUzMzAwNzggbCA0Ni44OTAwNjA3LDQ0LjM3OTQ5MjIgbCA0Ni41MDYyNzE2LDQzLjkxODk0NTMgbCA0Ni4wNDU3MjQ3LDQzLjUzNTE1NjMgbCA0NS42NjE5MzU3LDQzLjA3NDYwOTQgbCA0NS4yMDEzODg4LDQyLjY5MDgyMDMgbCA0NS4wNDc4NzMyLDQyLjA3Njc1NzggbCA0NS4wNDc4NzMyLDQxLjIzMjQyMTkgbCA0NS4yMDEzODg4LDQwLjU0MTYwMTYgbCA0NS42NjE5MzU3LDQwLjE1NzgxMjUgbCA0Ni4wNDU3MjQ3LDM5LjY5NzI2NTYgbCA0Ni41MDYyNzE2LDM5LjMxMzQ3NjYgbCA0Ni44OTAwNjA3LDM4Ljg1MjkyOTcgbCA0Ny41MDQxMjMyLDM4LjYyMjY1NjMgbCA0OC4xOTQ5NDM1LDM4LjQ2OTE0MDYgbCA0OC41Nzg3MzI2LDM4LjAwODU5MzggbCA0OS4wMzkyNzk0LDM3LjYyNDgwNDcgbCA0OS40MjMwNjg1LDM3LjE2NDI1NzggbCA0OS44ODM2MTU0LDM2Ljc4MDQ2ODggbCA1MC4yNjc0MDQ0LDM2LjMxOTkyMTkgbCA1MC43Mjc5NTEzLDM1LjkzNjEzMjggbCA1MC45NTgyMjQ3LDM1LjMyMjA3MDMgbCA1MC45NTgyMjQ3LDM0LjY5OTg2NjMgbCA1MC45NTgyMjQ3LDMzLjkxNTA0NDEgbCA1MC45NTgyMjQ3LDMyLjk0MDAxMDkgbCA1MC45NTgyMjQ3LDMyLjMwNzU1NDQgbCA1MC45NTgyMjQ3LDMxLjQ2ODMzNTIgbCA1MC45NTgyMjQ3LDMwLjU4OTEyMzIgbCA1MC45NTgyMjQ3LDI5LjgxMzg3NDIgbCA1MC45NTgyMjQ3LDI5LjAyMzY4NCBsIDUwLjk1ODIyNDcsMjguMTIzODA0NyBsIDUwLjk1ODIyNDcsMjcuMDAzOTIxMSBsIDUwLjk1ODIyNDcsMjYuMDM0Mzc1IGMgNTAuOTU4MjI0OSwyNS44NDYwMTI1IDUxLjExMTc0MDQsMjUuMzQzNTU0NyA1MS4xMTE3NDA0LDI1LjM0MzU1NDcgbCA1MS41NzIyODcyLDI0Ljk1OTc2NTYgbCA1MS44MDI1NjA3LDI0LjM0NTcwMzEgbCA1MS44MDI1NjA3LDIzLjUwMTM2NzIgbCA1MS41NzIyODcyLDIyLjgxMDU0NjkgbCA1MS4xMTE3NDA0LDIyLjQyNjc1NzggbCA1MC45NTgyMjQ3LDIxLjgxMjY5NTMgbCA1MC45NTgyMjQ3LDIxLjAyOTA5MTIgbCA1MC45NTgyMjQ3LDIwLjMzNTg4NTQgbCA1MC45NTgyMjQ3LDE5LjYzNDcxNjcgbCA1MC45NTgyMjQ3LDE4Ljk1NzUyNTggbCA1MC45NTgyMjQ3LDE4LjAwMDU2NTMgbCA1MC45NTgyMjQ3LDE3LjM1NTMxNDggbCA1MC45NTgyMjQ3LDE2Ljc0NjY3OTcgbCA1MC43Mjc5NTEzLDE2LjA1NTg1OTQgbCA1MC4yNjc0MDQ0LDE1LjY3MjA3MDMgbCA1MC4xMTM4ODg4LDE1LjA1ODAwNzggbCA0OS44ODM2MTU0LDE0LjM2NzE4NzUgbCA0OS40MjMwNjg1LDEzLjk4MzM5ODQgbCA0OS4wMzkyNzk0LDEzLjUyMjg1MTYgbCA0OC41Nzg3MzI2LDEzLjEzOTA2MjUgbCA0OC4xOTQ5NDM1LDEyLjY3ODUxNTYgbCA0Ny43MzQzOTY2LDEyLjI5NDcyNjYgbCA0Ny4zNTA2MDc2LDExLjgzNDE3OTcgbCA0Ni44OTAwNjA3LDExLjQ1MDM5MDYgbCA0Ni41MDYyNzE2LDEwLjk4OTg0MzggbCA0NS44MTU0NTEzLDEwLjc1OTU3MDMgbCA0NS4yMDEzODg4LDEwLjYwNjA1NDcgbCA0NS4yMDEzODg4LDEwLjE0NTUwNzggbCA0NS42NjE5MzU3LDkuNzYxNzE4NzUgbCA0NS44OTIyMDkxLDkuMTQ3NjU2MjUgbCA0NS44OTIyMDkxLDguMzAzMzIwMzEgbCA0Ni4wNDU3MjQ3LDcuNjEyNSBsIDQ2LjUwNjI3MTYsNy4yMjg3MTA5NCBsIDQ2Ljg5MDA2MDcsNi43NjgxNjQwNiBsIDQ3LjM1MDYwNzYsNi4zODQzNzUgbCA0Ny43MzQzOTY2LDUuOTIzODI4MTMgbCA0OC4zNDg0NTkxLDUuNjkzNTU0NjkgbCA0OS4wMzkyNzk0LDUuNTQwMDM5MDYgbCA0OS40MjMwNjg1LDUuMDc5NDkyMTkgbCA1MC4wMzcxMzEsNC44NDkyMTg3NSBsIDUwLjcyNzk1MTMsNS4wNzk0OTIxOSBsIDUxLjExMTc0MDQsNS41NDAwMzkwNiBsIDUxLjcyNTgwMjksNS42OTM1NTQ2OSBsIDUyLjU3MDEzODgsNS42OTM1NTQ2OSBsIDUzLjI2MDk1OTEsNS45MjM4MjgxMyBsIDUzLjY0NDc0ODIsNi4zODQzNzUgbCA1NC4yNTg4MTA3LDYuNTM3ODkwNjMgbCA1NC45NDk2MzEsNi43NjgxNjQwNiBsIDU1LjMzMzQyMDEsNy4yMjg3MTA5NCBsIDU1Ljc5Mzk2NjksNy42MTI1IGwgNTYuMTc3NzU2LDguMDczMDQ2ODggbCA1Ni42MzgzMDI5LDguNDU2ODM1OTQgbCA1Ny4wMjIwOTE5LDguOTE3MzgyODEgbCA1Ny40ODI2Mzg4LDkuMzAxMTcxODggbCA1Ny44NjY0Mjc5LDkuNzYxNzE4NzUgbCA1OC4zMjY5NzQ3LDEwLjE0NTUwNzggbCA1OC43MTA3NjM4LDEwLjYwNjA1NDcgbCA1OS4xNzEzMTA3LDEwLjYwNjA1NDcgbCA1OS41NTUwOTk3LDEwLjE0NTUwNzggbCA2MC4xNjkxNjIyLDkuOTE1MjM0MzggbCA2MS4wMTM0OTgyLDkuOTE1MjM0MzggbCA2MS43MDQzMTg1LDkuNzYxNzE4NzUgbCA2Mi4wODgxMDc2LDkuMzAxMTcxODggbCA2Mi41NDg2NTQ0LDguOTE3MzgyODEgbCA2Mi45MzI0NDM1LDguNDU2ODM1OTQgbCA2My4zOTI5OTA0LDguMDczMDQ2ODggbCA2My43NzY3Nzk0LDcuNjEyNSBsIDY0LjM5MDg0MTksNy4zODIyMjY1NiBsIDY1LjA4MTY2MjIsNy4yMjg3MTA5NCBsIDY1LjQ2NTQ1MTMsNi43NjgxNjQwNiBsIDY2LjA3OTUxMzgsNi41Mzc4OTA2MyBsIDY2LjY0NTY2NzksNi41Mzc4OTA2MyBsIDY3LjMwMjEwMiw2LjUzNzg5MDYzIGwgNjguMDU5MTg4Nyw2LjUzNzg5MDYzIGwgNjguNjEyNTIxNiw2LjUzNzg5MDYzIGwgNjkuMzAzMzQxOSw2LjM4NDM3NSBsIDY5LjY4NzEzMSw1LjkyMzgyODEzIGwgNzAuMzAxMTkzNSw1LjY5MzU1NDY5IGwgNzAuOTkwODE5Nyw1LjY5MzU1NDY5IGwgNzEuNTQ4MTIyMiw1LjY5MzU1NDY5IGwgNzIuMDQzNjAxNSw1LjY5MzU1NDY5IGwgNzIuNzg1ODM2NCw1LjY5MzU1NDY5IGwgNzMuNjc4NTM3Miw1LjY5MzU1NDY5IGwgNzQuMzY5MzU3Niw1LjkyMzgyODEzIGwgNzQuNzUzMTQ2Niw2LjM4NDM3NSBsIDc1LjM2NzIwOTEsNi41Mzc4OTA2MyBsIDc2LjA1ODAyOTQsNi43NjgxNjQwNiBsIDc2LjQ0MTgxODUsNy4yMjg3MTA5NCBsIDc3LjA1NTg4MSw3LjM4MjIyNjU2IGwgNzcuNzQ2NzAxMyw3LjYxMjUgbCA3OC4xMzA0OTA0LDguMDczMDQ2ODggbCA3OC41OTEwMzcyLDguNDU2ODM1OTQgbCA3OC45NzQ4MjYzLDguOTE3MzgyODEgbCA3OS40MzUzNzMyLDkuMzAxMTcxODggbCA3OS44MTkxNjIyLDkuNzYxNzE4NzUgbCA4MC4yNzk3MDkxLDEwLjE0NTUwNzggbCA4MC41MDk5ODI2LDEwLjgzNjMyODEgbCA4MC42NjM0OTgyLDExLjQ1MDM5MDYgbCA4MS4xMjQwNDUxLDExLjgzNDE3OTcgbCA4MS4zNTQzMTg1LDEyLjUyNSBsIDgxLjUwNzgzNDEsMTMuMTM5MDYyNSBsIDgxLjk2ODM4MSwxMy41MjI4NTE2IGwgODIuMTk4NjU0NCwxNC4yMTM2NzE5IGwgODIuMTk4NjU0NCwxNC45NDAwOTg5IGwgODIuMTk4NjU0NCwxNS45NzkxMDI1IGwgODIuMTk4NjU0NCwxNi43NDY2Nzk3IGwgODEuOTY4MzgxLDE3LjM2MDc0MjIgbCA4MS41MDc4MzQxLDE3Ljc0NDUzMTMgbCA4MS41MDc4MzQxLDE4LjIwNTA3ODEgbCA4MS45NjgzODEsMTguNTg4ODY3MiBsIDgxLjk2ODM4MSwxOS4wNDk0MTQxIGwgODEuNTA3ODM0MSwxOS40MzMyMDMxIGwgODEuMzU0MzE4NSwyMC4xMjQwMjM0IGwgODEuMzU0MzE4NSwyMC45NjUzODkzIGwgODEuMzU0MzE4NSwyMS44MTI2OTUzIGMgODEuMzU0MzE5MywyMi40MjY3NTgyIDgxLjUwNzgzNDEsMjIuNDI2NzU3OCA4MS41MDc4MzQxLDIyLjQyNjc1NzggbCA4MS45NjgzODEsMjIuODEwNTQ2OSBsIDgyLjE5ODY1NDQsMjMuNTAxMzY3MiBsIDgyLjE5ODY1NDQsMjQuMzQyMzA3NiBsIDgyLjE5ODY1NDQsMjUuMjE5MjgyOCBsIDgyLjE5ODY1NDQsMjUuOTI5NTc3MiBsIDgyLjE5ODY1NDQsMjcuMDcxMzgwNyBsIDgyLjE5ODY1NDQsMjcuNjg0NjkwOCBsIDgyLjE5ODY1NDQsMjguNTY3MzgyOCBjIDgyLjE5ODY1MjksMjkuMTgxNDQ0NSA4Mi4zNTIxNzAxLDI5LjE4MTQ0NTMgODIuMzUyMTcwMSwyOS4xODE0NDUzIGwgODIuODEyNzE2OSwyOS41NjUyMzQ0IGwgODMuMDQyOTkwNCwzMC4yNTYwNTQ3IGwgODMuMDQyOTkwNCwzMS4xMDAzOTA2IGwgODIuODEyNzE2OSwzMS43MTQ0NTMxIGwgODIuMzUyMTcwMSwzMi4wOTgyNDIyIGwgODIuMTk4NjU0NCwzMi43ODkwNjI1IGwgODIuMzUyMTcwMSwzMy40MDMxMjUgbCA4Mi44MTI3MTY5LDMzLjc4NjkxNDEgbCA4My4wNDI5OTA0LDM0LjQ3NzczNDQgbCA4My4xOTY1MDYsMzUuMDkxNzk2OSBsIDgzLjY1NzA1MjksMzUuNDc1NTg1OSBsIDg0LjA0MDg0MTksMzUuOTM2MTMyOCBsIDg0LjY1NDkwNDQsMzYuMDg5NjQ4NCBsIDg1LjM0NTcyNDcsMzYuMzE5OTIxOSBsIDg1LjcyOTUxMzgsMzYuNzgwNDY4OCBsIDg2LjM0MzU3NjMsMzYuOTMzOTg0NCBsIDg3LjAzNDM5NjYsMzcuMTY0MjU3OCBsIDg3LjI2NDY3MDEsMzcuODU1MDc4MSBsIDg3LjQxODE4NTcsMzguNDY5MTQwNiBsIDg3Ljg3ODczMjYsMzguODUyOTI5NyBsIDg4LjEwOTAwNiwzOS41NDM3NSBsIDg3Ljg3ODczMjYsNDAuMTU3ODEyNSBsIDg3LjQxODE4NTcsNDAuNTQxNjAxNiB6IG0gMTMwLjQ4ODkyNCwzMi4wOTgyNDIyIGwgMTMwLjMzNTQwOCwzMi43ODkwNjI1IGwgMTMwLjMzNTQwOCwzMy42ODg2NjA2IGwgMTMwLjMzNTQwOCwzNC40Nzc3MzQ0IGwgMTMwLjEwNTEzNSwzNS4wOTE3OTY5IGwgMTI5LjY0NDU4OCwzNS40NzU1ODU5IGwgMTI5LjQ5MTA3MiwzNi4xNjY0MDYzIGwgMTI5LjQ5MTA3MiwzNi45NzQyNjMxIGwgMTI5LjQ5MTA3MiwzNy44NTUwNzgxIGwgMTI5LjI2MDc5OSwzOC40NjkxNDA2IGwgMTI4LjgwMDI1MiwzOC44NTI5Mjk3IGwgMTI4LjY0NjczNiwzOS41NDM3NSBsIDEyOC40MTY0NjMsNDAuMTU3ODEyNSBsIDEyNy45NTU5MTYsNDAuNTQxNjAxNiBsIDEyNy41NzIxMjcsNDEuMDAyMTQ4NCBsIDEyNy4xMTE1OCw0MS4zODU5Mzc1IGwgMTI2LjcyNzc5MSw0MS44NDY0ODQ0IGwgMTI2LjI2NzI0NCw0Mi4yMzAyNzM0IGwgMTI1Ljg4MzQ1NSw0Mi42OTA4MjAzIGwgMTI1LjQyMjkwOCw0My4wNzQ2MDk0IGwgMTI1LjAzOTExOSw0My41MzUxNTYzIGwgMTI0LjU3ODU3Miw0My45MTg5NDUzIGwgMTI0LjE5NDc4Myw0NC4zNzk0OTIyIGwgMTIzLjUwMzk2Myw0NC41MzMwMDc4IGwgMTIyLjg4OTksNDQuNzYzMjgxMyBsIDEyMi41MDYxMTEsNDUuMjIzODI4MSBsIDEyMS44MTUyOTEsNDUuMzc3MzQzOCBsIDEyMC45NzA5NTUsNDUuMzc3MzQzOCBsIDEyMC4zNTY4OTIsNDUuNjA3NjE3MiBsIDExOS45NzMxMDMsNDYuMDY4MTY0MSBsIDExOS4yODIyODMsNDYuMjIxNjc5NyBsIDExOC42MjQxODUsNDYuMjIxNjc5NyBsIDExNy45NDE1NjYsNDYuMjIxNjc5NyBsIDExNy4xMTg1MDcsNDYuMjIxNjc5NyBsIDExNi40NDIwNDMsNDYuMjIxNjc5NyBsIDExNS43MDc1NDgsNDYuMjIxNjc5NyBsIDExNS4wNjA2MDMsNDYuMjIxNjc5NyBsIDExNC40NDY1NDEsNDYuMDY4MTY0MSBsIDExNC4wNjI3NTIsNDUuNjA3NjE3MiBsIDExMy4zNzE5MzIsNDUuMzc3MzQzOCBsIDExMi41Mjc1OTYsNDUuMzc3MzQzOCBsIDExMS45MTM1MzMsNDUuMjIzODI4MSBsIDExMS41Mjk3NDQsNDQuNzYzMjgxMyBsIDExMC44Mzg5MjQsNDQuNTMzMDA3OCBsIDExMC4yMjQ4NjEsNDQuMzc5NDkyMiBsIDEwOS44NDEwNzIsNDMuOTE4OTQ1MyBsIDEwOS4zODA1MjUsNDMuNTM1MTU2MyBsIDEwOC45OTY3MzYsNDMuMDc0NjA5NCBsIDEwOC41MzYxODksNDIuNjkwODIwMyBsIDEwOC4xNTI0LDQyLjIzMDI3MzQgbCAxMDcuNjkxODUzLDQxLjg0NjQ4NDQgbCAxMDcuMzA4MDY0LDQxLjM4NTkzNzUgbCAxMDYuODQ3NTE3LDQxLjAwMjE0ODQgbCAxMDYuNDYzNzI4LDQwLjU0MTYwMTYgbCAxMDYuMDAzMTgyLDQwLjE1NzgxMjUgbCAxMDUuNjE5MzkyLDM5LjY5NzI2NTYgbCAxMDUuMTU4ODQ2LDM5LjMxMzQ3NjYgbCAxMDUuMDA1MzMsMzguNjk5NDE0MSBsIDEwNC43NzUwNTcsMzguMDA4NTkzOCBsIDEwNC4zMTQ1MSwzNy42MjQ4MDQ3IGwgMTA0LjE2MDk5NCwzNy4wMTA3NDIyIGwgMTA0LjE2MDk5NCwzNi4xNjY0MDYzIGwgMTAzLjkzMDcyMSwzNS40NzU1ODU5IGwgMTAzLjQ3MDE3NCwzNS4wOTE3OTY5IGwgMTAzLjMxNjY1OCwzNC40Nzc3MzQ0IGwgMTAzLjMxNjY1OCwzMy43NDQwNTM0IGwgMTAzLjMxNjY1OCwzMi43ODkwNjI1IGMgMTAzLjMxNjY1OCwzMi40MzExNzA3IDEwMy40NzAxNzQsMzIuMDk4MjQyMiAxMDMuNDcwMTc0LDMyLjA5ODI0MjIgbCAxMDMuOTMwNzIxLDMxLjcxNDQ1MzEgbCAxMDQuMTYwOTk0LDMxLjEwMDM5MDYgbCAxMDQuMTYwOTk0LDMwLjM5NjI2MzkgbCAxMDQuMTYwOTk0LDI5LjUxNTA5NDEgbCAxMDQuMTYwOTk0LDI4LjQzMDY1NzQgbCAxMDQuMTYwOTk0LDI3LjcyMzA0NjkgYyAxMDQuMTYwOTk0LDI3LjU0NDI5MTcgMTA0LjMxNDUxLDI3LjAzMjIyNjYgMTA0LjMxNDUxLDI3LjAzMjIyNjYgbCAxMDQuNzc1MDU3LDI2LjY0ODQzNzUgbCAxMDUuMDA1MzMsMjYuMDM0Mzc1IGwgMTA0Ljc3NTA1NywyNS4zNDM1NTQ3IGwgMTA0LjMxNDUxLDI0Ljk1OTc2NTYgbCAxMDQuMTYwOTk0LDI0LjM0NTcwMzEgbCAxMDQuMzE0NTEsMjMuNjU0ODgyOCBsIDEwNC43NzUwNTcsMjMuMjcxMDkzOCBsIDEwNS4wMDUzMywyMi42NTcwMzEzIGwgMTA1LjAwNTMzLDIxLjU0MTE2MzQgbCAxMDUuMDA1MzMsMjAuNjI3MTQxNiBsIDEwNS4wMDUzMywxOS45OTIyNDIyIGwgMTA1LjAwNTMzLDE5LjI1NjcwODUgbCAxMDUuMDA1MzMsMTguNDk1NTk2NSBsIDEwNS4wMDUzMywxNy42NDc3NTg0IGwgMTA1LjAwNTMzLDE2LjgwMDAwMDIgbCAxMDUuMDA1MzMsMTUuOTAyMzQzOCBsIDEwNC43NzUwNTcsMTUuMjExNTIzNCBsIDEwNC4zMTQ1MSwxNC44Mjc3MzQ0IGwgMTA0LjE2MDk5NCwxNC4yMTM2NzE5IGwgMTAzLjkzMDcyMSwxMy41MjI4NTE2IGwgMTAzLjQ3MDE3NCwxMy4xMzkwNjI1IGwgMTAzLjA4NjM4NSwxMi42Nzg1MTU2IGwgMTAyLjM5NTU2NCwxMi40NDgyNDIyIGwgMTAxLjc4MTUwMiwxMi4yOTQ3MjY2IGwgMTAxLjM5NzcxMywxMS44MzQxNzk3IGwgMTAwLjcwNjg5MiwxMS42MDM5MDYzIGwgMTAwLjA1MzM1NywxMS42MDM5MDYzIGwgOTkuNDI4ODQ4MywxMS42MDM5MDYzIGwgOTguODIwMDA2NywxMS42MDM5MDYzIGwgOTguMTczODg0NiwxMS42MDM5MDYzIGwgOTcuNTU5ODIyMSwxMS40NTAzOTA2IGwgOTcuMTc2MDMzMSwxMC45ODk4NDM4IGwgOTYuNzE1NDg2MiwxMC42MDYwNTQ3IGwgOTYuMzMxNjk3MSwxMC4xNDU1MDc4IGwgOTUuODcxMTUwMyw5Ljc2MTcxODc1IGwgOTUuNDg3MzYxMiw5LjMwMTE3MTg4IGwgOTUuMDI2ODE0Myw4LjkxNzM4MjgxIGwgOTQuODczMjk4Nyw4LjMwMzMyMDMxIGwgOTUuMDI2ODE0Myw3LjYxMjUgbCA5NS40ODczNjEyLDcuMjI4NzEwOTQgbCA5NS44NzExNTAzLDYuNzY4MTY0MDYgbCA5Ni4zMzE2OTcxLDYuMzg0Mzc1IGwgOTYuNzE1NDg2Miw1LjkyMzgyODEzIGwgOTcuMTc2MDMzMSw1LjU0MDAzOTA2IGwgOTcuNTU5ODIyMSw1LjA3OTQ5MjE5IGwgOTguMTczODg0Niw0Ljg0OTIxODc1IGwgOTguNzMzNzYwMSw0Ljg0OTIxODc1IGwgOTkuMjE4Mzg3NSw0Ljg0OTIxODc1IGwgOTkuNzAxOTc1OCw0Ljg0OTIxODc1IGwgMTAwLjE3ODUzLDQuODQ5MjE4NzUgbCAxMDAuNjA0MDg4LDQuODQ5MjE4NzUgbCAxMDEuNDAzODA3LDQuODQ5MjE4NzUgbCAxMDIuMDAyNjU3LDQuODQ5MjE4NzUgbCAxMDIuMzk1NTY0LDQuODQ5MjE4NzUgbCAxMDMuMDg2Mzg1LDUuMDc5NDkyMTkgbCAxMDMuNDcwMTc0LDUuNTQwMDM5MDYgbCAxMDMuOTMwNzIxLDUuNTQwMDM5MDYgbCAxMDQuMzE0NTEsNS4wNzk0OTIxOSBsIDEwNC43NzUwNTcsNC42OTU3MDMxMyBsIDEwNS4wMDUzMyw0LjA4MTY0MDYzIGwgMTA1LjE1ODg0NiwzLjM5MDgyMDMxIGwgMTA1LjYxOTM5MiwzLjAwNzAzMTI1IGwgMTA1Ljg0OTY2NiwyLjM5Mjk2ODc1IGwgMTA2LjAwMzE4MiwxLjcwMjE0ODQ0IGwgMTA2LjYxNzI0NCwxLjQ3MTg3NSBsIDEwNy4zMDgwNjQsMS4zMTgzNTkzOCBsIDEwNy41MzgzMzgsMC43MDQyOTY4NzUgbCAxMDcuNjkxODUzLDAuMDEzNDc2NTYyNSBsIDEwOC4xNTI0LDAuMDEzNDc2NTYyNSBsIDEwOC4zODI2NzQsMC43MDQyOTY4NzUgbCAxMDguNTM2MTg5LDEuMzE4MzU5MzggbCAxMDguOTk2NzM2LDEuNzAyMTQ4NDQgbCAxMDkuMzgwNTI1LDIuMTYyNjk1MzEgbCAxMDkuODQxMDcyLDIuNTQ2NDg0MzggbCAxMTAuMDcxMzQ2LDMuMjM3MzA0NjkgbCAxMTAuMDcxMzQ2LDQuMDgxNjQwNjMgbCAxMTAuMjI0ODYxLDQuNjk1NzAzMTMgbCAxMTAuNjg1NDA4LDUuMDc5NDkyMTkgbCAxMTEuMDY5MTk3LDUuNTQwMDM5MDYgbCAxMTEuNjgzMjYsNS42OTM1NTQ2OSBsIDExMi4zNzQwOCw1LjkyMzgyODEzIGwgMTEyLjc1Nzg2OSw2LjM4NDM3NSBsIDExMy4yMTg0MTYsNi4zODQzNzUgbCAxMTMuNjAyMjA1LDUuOTIzODI4MTMgbCAxMTQuMjE2MjY3LDUuNjkzNTU0NjkgbCAxMTUuMDA0NjI3LDUuNjkzNTU0NjkgbCAxMTUuOTA0OTM5LDUuNjkzNTU0NjkgbCAxMTYuNTk1NzYsNS45MjM4MjgxMyBsIDExNi45Nzk1NDksNi4zODQzNzUgbCAxMTcuNTkzNjExLDYuNTM3ODkwNjMgbCAxMTguMjUxNzgyLDYuNTM3ODkwNjMgbCAxMTguNzMyMDkzLDYuNTM3ODkwNjMgbCAxMTkuMzg5MDU0LDYuNTM3ODkwNjMgYyAxMTkuNjEwMDkyLDYuNTM3ODkwNjMgMTE5LjgzMTEzLDYuNTM3ODkwNjMgMTIwLjA1MjE2OSw2LjUzNzg5MDYzIGMgMTIwLjIwMzc2NSw2LjUzNzg5MDYzIDEyMC4zODE4MDMsNi41Mzc4OTA2MyAxMjAuNTcyNDQ2LDYuNTM3ODkwNjMgYyAxMjAuNzY0OTIyLDYuNTM3ODkwNjMgMTIwLjk1NzM5OCw2LjUzNzg5MDYzIDEyMS4xNDk4NzQsNi41Mzc4OTA2MyBjIDEyMS4zNTIxMzUsNi41Mzc4OTA2MyAxMjEuNjM1Njg2LDYuNTM3ODkwNjMgMTIxLjkwMjUxMyw2LjUzNzg5MDYzIGwgMTIyLjY1OTYyNyw2LjUzNzg5MDYzIGMgMTIyLjg5NzU3LDYuNTM3ODkwNjMgMTIzLjM1MDQ0Nyw2Ljc2ODE2NDA2IDEyMy4zNTA0NDcsNi43NjgxNjQwNiBsIDEyMy43MzQyMzYsNy4yMjg3MTA5NCBsIDEyNC4xOTQ3ODMsNy42MTI1IGwgMTI0LjU3ODU3Miw4LjA3MzA0Njg4IGwgMTI1LjAzOTExOSw4LjQ1NjgzNTk0IGwgMTI1LjI2OTM5Miw5LjE0NzY1NjI1IGwgMTI1LjAzOTExOSw5Ljc2MTcxODc1IGwgMTI0LjU3ODU3MiwxMC4xNDU1MDc4IGwgMTI0LjE5NDc4MywxMC42MDYwNTQ3IGwgMTIzLjczNDIzNiwxMC45ODk4NDM4IGwgMTIzLjM1MDQ0NywxMS40NTAzOTA2IGwgMTIyLjY1OTYyNywxMS42MDM5MDYzIGwgMTIxLjk4ODIsMTEuNjAzOTA2MyBsIDEyMS4zNDkzODQsMTEuNjAzOTA2MyBsIDEyMC43ODA5ODgsMTEuNjAzOTA2MyBjIDEyMC41OTM2MDEsMTEuNjAzOTA2MyAxMjAuNDA2MjE0LDExLjYwMzkwNjMgMTIwLjIxODgyNywxMS42MDM5MDYzIGwgMTE5LjQ0MjQ0OCwxMS42MDM5MDYzIGwgMTE4Ljc3NDEzNywxMS42MDM5MDYzIGwgMTE4LjE3NjI0NiwxMS42MDM5MDYzIGwgMTE3LjUyMDI0NSwxMS42MDM5MDYzIGwgMTE2LjgzNjY2NywxMS42MDM5MDYzIGwgMTE2LjIzODc3NiwxMS42MDM5MDYzIGMgMTE2LjA1MTg2MywxMS42MDM5MDYzIDExNS44Njg2NDksMTEuNjAzOTA2MyAxMTUuNjkwOTIyLDExLjYwMzkwNjMgYyAxMTUuNDkwNjIyLDExLjYwMzkwNjMgMTE1LjI5NzI5MiwxMS42MDM5MDYzIDExNS4xMTM0OTQsMTEuNjAzOTA2MyBsIDExNC41MzYwNjYsMTEuNjAzOTA2MyBjIDExNC4zNzQ5MDcsMTEuNjAzOTA2MyAxMTQuMjI2NDQxLDExLjYwMzkwNjMgMTE0LjA5MzMyMiwxMS42MDM5MDYzIGMgMTEzLjY0NTk2NywxMS42MDM5MDYzIDExMy4zNzE5MzIsMTEuNjAzOTA2MyAxMTMuMzcxOTMyLDExLjYwMzkwNjMgbCAxMTIuNzU3ODY5LDExLjgzNDE3OTcgbCAxMTIuMzc0MDgsMTIuMjk0NzI2NiBsIDExMS45MTM1MzMsMTIuNjc4NTE1NiBsIDExMS41Mjk3NDQsMTMuMTM5MDYyNSBsIDExMS4wNjkxOTcsMTMuNTIyODUxNiBsIDExMC42ODU0MDgsMTMuOTgzMzk4NCBsIDExMC4yMjQ4NjEsMTQuMzY3MTg3NSBsIDExMC4wNzEzNDYsMTUuMDU4MDA3OCBsIDExMC4wNzEzNDYsMTUuOTkwMjg5NiBsIDExMC4wNzEzNDYsMTYuNzQ2Njc5NyBjIDExMC4wNzEzNDUsMTcuMjE4OTIyOCAxMTAuMjI0ODYxLDE3LjM2MDc0MjIgMTEwLjIyNDg2MSwxNy4zNjA3NDIyIGwgMTEwLjY4NTQwOCwxNy43NDQ1MzEzIGwgMTEwLjY4NTQwOCwxOC4yMDUwNzgxIGwgMTEwLjIyNDg2MSwxOC41ODg4NjcyIGwgMTEwLjA3MTM0NiwxOS4yNzk2ODc1IGwgMTEwLjA3MTM0NiwyMC4xNzMwNDgzIGMgMTEwLjA3MTM0NiwyMC4zOTUwNTUyIDExMC4wNzEzNDYsMjAuNjM0Njc2NiAxMTAuMDcxMzQ2LDIwLjg4NzcxOTggbCAxMTAuMDcxMzQ2LDIxLjY1Njk4NDggbCAxMTAuMDcxMzQ2LDIyLjM5MDUyMDIgbCAxMTAuMDcxMzQ2LDIzLjE4NzM2MTggbCAxMTAuMDcxMzQ2LDIzLjk3MDkzNDYgbCAxMTAuMDcxMzQ2LDI0LjY1MzU1MzQgbCAxMTAuMDcxMzQ2LDI1LjQ3MzczNTEgbCAxMTAuMDcxMzQ2LDI2LjI2NTM4MSBsIDExMC4wNzEzNDYsMjYuOTEzMjI5NCBsIDExMC4wNzEzNDYsMjcuNTQ0MjkyMSBjIDExMC4wNzEzNDYsMjguMTc0MjU0MiAxMTAuMDcxMzQ2LDI4LjU2NzM4MjggMTEwLjA3MTM0NiwyOC41NjczODI4IGwgMTEwLjIyNDg2MSwyOS4xODE0NDUzIGwgMTEwLjY4NTQwOCwyOS41NjUyMzQ0IGwgMTEwLjkxNTY4MiwzMC4yNTYwNTQ3IGwgMTEwLjkxNTY4MiwzMS4xMDAzOTA2IGwgMTExLjA2OTE5NywzMS43MTQ0NTMxIGwgMTExLjUyOTc0NCwzMi4wOTgyNDIyIGwgMTExLjc2MDAxNywzMi43ODkwNjI1IGwgMTExLjc2MDAxNywzMy42MzMzOTg0IGwgMTExLjkxMzUzMywzNC4yNDc0NjA5IGwgMTEyLjM3NDA4LDM0LjYzMTI1IGwgMTEyLjYwNDM1MywzNS4zMjIwNzAzIGwgMTEyLjc1Nzg2OSwzNS45MzYxMzI4IGwgMTEzLjM3MTkzMiwzNi4wODk2NDg0IGwgMTE0LjA2Mjc1MiwzNi4zMTk5MjE5IGwgMTE0LjQ0NjU0MSwzNi43ODA0Njg4IGwgMTE1LjA2MDYwMywzNi45MzM5ODQ0IGwgMTE1Ljc0MjMxOCwzNi45MzM5ODQ0IGwgMTE2LjIyNjg2NiwzNi45MzM5ODQ0IGwgMTE2LjY4MzkxNywzNi45MzM5ODQ0IGwgMTE3LjA0MTg1MiwzNi45MzM5ODQ0IGwgMTE3LjcyNzM0OCwzNi45MzM5ODQ0IGwgMTE4LjQzNzk0NywzNi45MzM5ODQ0IGwgMTE5LjEyODc2NywzNi43ODA0Njg4IGwgMTE5LjUxMjU1NywzNi4zMTk5MjE5IGwgMTE5Ljk3MzEwMywzNS45MzYxMzI4IGwgMTIwLjM1Njg5MiwzNS40NzU1ODU5IGwgMTIwLjgxNzQzOSwzNS4wOTE3OTY5IGwgMTIxLjA0NzcxMywzNC40Nzc3MzQ0IGwgMTIxLjIwMTIyOCwzMy43ODY5MTQxIGwgMTIxLjY2MTc3NSwzMy40MDMxMjUgbCAxMjEuODkyMDQ5LDMyLjc4OTA2MjUgbCAxMjEuODkyMDQ5LDMxLjk0NDcyNjYgbCAxMjIuMDQ1NTY0LDMxLjI1MzkwNjMgbCAxMjIuNTA2MTExLDMwLjg3MDExNzIgbCAxMjIuNzM2Mzg1LDMwLjI1NjA1NDcgbCAxMjIuODg5OSwyOS41NjUyMzQ0IGwgMTIzLjM1MDQ0NywyOS4xODE0NDUzIGwgMTIzLjczNDIzNiwyOC43MjA4OTg0IGwgMTI0LjM0ODI5OSwyOC40OTA2MjUgbCAxMjUuMDM5MTE5LDI4LjMzNzEwOTQgbCAxMjUuNDIyOTA4LDI3Ljg3NjU2MjUgbCAxMjYuMDM2OTcxLDI3LjY0NjI4OTEgbCAxMjYuODgxMzA3LDI3LjY0NjI4OTEgbCAxMjcuNTcyMTI3LDI3Ljg3NjU2MjUgbCAxMjcuOTU1OTE2LDI4LjMzNzEwOTQgbCAxMjguNTY5OTc4LDI4LjQ5MDYyNSBsIDEyOS4yNjA3OTksMjguNzIwODk4NCBsIDEyOS42NDQ1ODgsMjkuMTgxNDQ1MyBsIDEzMC4xMDUxMzUsMjkuNTY1MjM0NCBsIDEzMC4zMzU0MDgsMzAuMjU2MDU0NyBsIDEzMC40ODg5MjQsMzAuODcwMTE3MiBsIDEzMC45NDk0NzEsMzEuMjUzOTA2MyBsIDEzMC45NDk0NzEsMzEuNzE0NDUzMSBsIDEzMC40ODg5MjQsMzIuMDk4MjQyMiB6IFwiO1xudmFyIGxldHRyZVMgPSBcIm0gMTM4LjgxMTUxNywzNS4xMDYxMDExIGMgMTM4LjYzNzgxLDM1LjczNjUxMTEgMTM4LjYwNjU2NCwzNS44MDEwNjE0IDEzOC4xOTg0NTUsMzYuMjA5MTcwNCBjIDEzOC4wNTY3NTMsMzYuMzUwODcyMiAxMzcuOTk5MTU1LDM3LjAyMTU5MjggMTM3Ljk0ODc5OCwzNy4yMjMwMTk5IGMgMTM3Ljg2OTUzNywzNy41NDAwNjIyIDEzNy45MjY0MTEsMzcuODI5MjMzMSAxMzcuOTQ4Nzk4LDM4LjA3NzU0NzUgYyAxMzcuOTgzMDc4LDM4LjQ1Nzc4MDYgMTM3LjkyNzc2NiwzOC45NDE1MzY3IDEzOC4xMjYwMzcsMzkuMjIxNzQ3OSBjIDEzOC4yNDA3MTMsMzkuMzgzODE2NCAxMzguNDM4NTk2LDM5LjQzMzgzNSAxMzguNTYwNTQ0LDM5LjYyMDE5NjQgYyAxMzguNzQ3MDY2LDM5LjkwNTI0MDggMTM4Ljc5NTkxLDQwLjQ1ODQ1MjcgMTM4Ljk1MTU5OSw0MC44MzM5MzMgYyAxMzkuMDYxMTMxLDQxLjA5ODA5MjIgMTM5LjI2MzkxOCw0MS4yMDQ1OTU4IDEzOS40Mjk1NTYsNDEuMzcwMjM0MyBjIDEzOS41OTI2NzIsNDEuNTMzMzQ5NSAxMzkuNjI2NDEyLDQxLjkzNDUzNjMgMTM5LjYyMjE2OSw0Mi40MjQ1ODc2IGMgMTM5LjYxODMyLDQyLjg2OTA2MjUgMTM5LjU4MzIyNiw0My4zODY2NDEgMTM5LjU4NTA3Niw0My44NjU5OTQ3IGMgMTM5LjU4NjU3OSw0NC4yNTU0NTk0IDEzOS42MTI0NzEsNDQuNjE5NjkxMiAxMzkuNjk5MzI0LDQ0Ljg5ODk4MDYgYyAxMzkuNzk4NTE4LDQ1LjIxNzk1NTcgMTQwLjA1ODUzMyw0NS4zMzY4NTE1IDE0MC4yNDkzODIsNDUuNTE2MjQzOCBjIDE0MC40NjgzNzIsNDUuNzIyMDg3NyAxNDAuNjY3Njk4LDQ1Ljk1Mjc0MDUgMTQwLjg4MzcxNCw0Ni4yMDY0MTk2IGMgMTQxLjI4Nzk4Nyw0Ni42ODExNzkxIDE0MS4zNjM1NjMsNDYuODA0MjY2MiAxNDEuNzU4NzMzLDQ2Ljg4ODk0NTIgYyAxNDIuMTUzOTAzLDQ2Ljk3MzYyNDEgMTQyLjM4OTU3MSw0Ni44ODg5NDUyIDE0Mi41NzczLDQ2Ljg4ODk0NTIgYyAxNDIuOTE5MDE4LDQ2Ljg4ODk0NTIgMTQzLjA2MjU4Niw0Ni41NDQ3OTI0IDE0My4xNzAwNTUsNDYuNDM3MzIyOSBjIDE0My41MjQ0LDQ2LjA4Mjk3NzYgMTQ0LjMwMjYxLDQ2LjA4NDkzMjcgMTQ0LjQ1NjEwNyw0NS45NTc2OTI4IGMgMTQ0LjYwOTYwMyw0NS44MzA0NTI4IDE0NC41MzgyNTIsNDUuNzg4ODk4IDE0NC42ODAxNjksNDUuNjQ2OTgxMSBjIDE0NC44NjczOTEsNDUuNDU5NzU5MiAxNDUuMzQwOCw0NS4zNDkyMzk1IDE0NS45MjEwNDcsNDUuMjkxOTMyOSBjIDE0Ni40MTY4NzksNDUuMjQyOTYzMyAxNDYuOTkwNzIzLDQ1LjIzMjg1IDE0Ny41MzA2NzMsNDUuMjQ2OTM2NSBjIDE0OC4yMzQ1OTksNDUuMjY1MzAwOCAxNDguODgwOTE3LDQ1LjMyNDc5NTIgMTQ5LjIyMTY2MSw0NS4zOTI5NDQyIGMgMTQ5Ljc5MzUyNyw0NS41MDczMTcyIDE0OS41OTIyOTMsNDUuOTQwNjc5NCAxNTAuMzExMzQ0LDQ2LjA4NDQ4OTYgYyAxNTAuNTkwNjE1LDQ2LjE0MDM0MzcgMTUwLjk4MzUwOCw0Ni4xMzI4NTk2IDE1MS4zODU3MDUsNDYuMTMzOTgzMyBjIDE1MS43NDA0NjEsNDYuMTM0OTc0NSAxNTIuMTAyNDU1LDQ2LjE0MjY2MjYgMTUyLjQwMDEwMSw0Ni4yMDY0MTk2IGMgMTUzLjAzNTE5NSw0Ni4zNDI0NTk2IDE1Mi44NDMyNTksNDYuNjk4OTY0OSAxNTMuMzAzNzk3LDQ2Ljg1NTM1NjMgYyAxNTMuNzQwMTA5LDQ3LjAwMzUyMDggMTU0LjEyNjY1OSw0Ny4wNjA3MTc3IDE1NC40NjI3NjcsNDcuMDUyNTMwNyBjIDE1NC44NDg4NTgsNDcuMDQzMTI2MyAxNTUuMTY4Mzg5LDQ2Ljk0NzQ0NTEgMTU1LjQyMDMyOSw0Ni44MDQyNjYyIGMgMTU1LjcwMTE1Myw0Ni42NDQ2NzIzIDE1NS44NTY1ODIsNDYuMzAxODg1MiAxNTYuMjM4NDQ0LDQ2LjIwNjQxOTYgYyAxNTYuNjM3NTUsNDYuMTA2NjQzIDE1Ny4wNjI0NjMsNDYuMDc3Mzg2NCAxNTcuNTAxMTkzLDQ2LjA3OTkwMDMgYyAxNTcuOTYyNTE0LDQ2LjA4MjU0MzYgMTU4LjQzOTExMSw0Ni4xMjAzMTMzIDE1OC45MTcwNDcsNDYuMTQ4MTYwOSBjIDE1OS4zOTUxNTgsNDYuMTc2MDE4NyAxNTkuODc0NjA5LDQ2LjE5Mzk0NyAxNjAuMzQxNDQ1LDQ2LjE1Njg0NzcgYyAxNjAuNjk1NjI0LDQ2LjEyODcwMTEgMTYxLjA0MjU0Myw0Ni4wNjg4ODA5IDE2MS4zNzYxMDcsNDUuOTU3NjkyOCBjIDE2MS41NTYxNjUsNDUuODk3NjczNCAxNjEuNTk1NTAxLDQ1LjY1NDIxNTkgMTYxLjg1NTk1Nyw0NS41MTYyNDM4IGMgMTYyLjIxNDM1Nyw0NS4zMjYzODcyIDE2Mi43NjA0NSw0NS4yMDQ1NzIyIDE2Mi45OTkxMjcsNDUuMTI1MDEzMyBjIDE2My4zMDE0NjksNDUuMDI0MjMyNiAxNjMuMjU0MjExLDQ0LjgxNzY4NDMgMTYzLjU0Mjc5OCw0NC42NzMzOTEgYyAxNjMuODExOTIxLDQ0LjUzODgyOTMgMTY0LjE5ODc1MSw0NC41MDQwMjk0IDE2NC41Mzc0NjksNDQuMzYyODk5MSBjIDE2NC44NzYxODcsNDQuMjIxNzY4NyAxNjUuMjQ1NjI4LDQzLjc5NTg3MTggMTY1LjU5NTk2LDQzLjQ0NTU0MDcgYyAxNjUuNjk0MDc0LDQzLjM0NzQyNjMgMTY2LjA4OTQyNyw0Mi45NTE1NzcgMTY2LjQxNDAzMSw0Mi41ODQ2MzM3IGMgMTY2LjczODYzNSw0Mi4yMTc2OTAzIDE2Ni42NDkzOTgsNDEuNjk1ODM4MSAxNjYuODIzODEsNDEuMzcwMjM0MyBjIDE2Ni45NDI1ODcsNDEuMTQ4NDkyOCAxNjcuMjMxNzI4LDQxLjA3MTcyMDIgMTY3LjM0NjE2NCw0MC44MzM5MzMgYyAxNjcuNTE0Nzk3LDQwLjQ4MzUyODYgMTY3LjQ3MDU2NCwzOS45NTg0MzM5IDE2Ny42NDI1NDIsMzkuNjIwMTk2NCBjIDE2Ny43NTIzMzIsMzkuNDA0MjY4MSAxNjguMDY0NzU2LDM5LjMzMDQ2NDYgMTY4LjE2NDczMiwzOS4xMjYyMzI4IGMgMTY4LjMzNzU3NiwzOC43NzMxNDMyIDE2OC4zODIyMzUsMzguMjE5MjQyNSAxNjguNTE3Mzk1LDM3Ljk1NDgzNzUgYyAxNjguNjYzODAyLDM3LjY2ODQyOTQgMTY4Ljk0NjM1MywzNy41MzY2OTUxIDE2OS4wMTEzNTgsMzcuNDA2Njg1NyBjIDE2OS4xNjE5OSwzNy4xMDU0MjE0IDE2OS4xOTY0MTgsMzYuNTgzNjE5NCAxNjkuMjE2MDAxLDM2LjA1MzQzMDEgYyAxNjkuMjM2ODMsMzUuNDg5NTA1NSAxNjkuMjQwODY2LDM0LjkxNjA5MjQgMTY5LjM1MDA3NCwzNC41ODg0NjkxIGMgMTY5LjQ2OTc3NCwzNC4yMjkzNjkzIDE2OS44MTU4MSwzNC4yOTA2OTAyIDE2OS45NzEwNTgsMzMuNzY4NTAyNyBjIDE3MC4xMjYzMDYsMzMuMjQ2MzE1MiAxNjkuNjc4NjEzLDMyLjc2NzYzNTQgMTY5LjQwOTI5LDMyLjQ5ODMxMjkgYyAxNjkuMjA5MzA1LDMyLjI5ODMyNzMgMTY5LjM0MTAzNSwzMS40NDM3NjUzIDE2OC41NzM4NSwzMC44NDcwNjU5IGMgMTY4LjMxOTgxMSwzMC42NDk0ODA1IDE2OC40NTY1NTUsMzAuMDA4NjY2NyAxNjguMTY0NzMyLDI5LjYwNTMyMzQgYyAxNjcuOTkzNDM3LDI5LjM2ODU2OTEgMTY3LjcyNDE2OSwyOS4yODU2NzM5IDE2Ny42NDI1NDIsMjkuMDQwNzk0NiBjIDE2Ny40NTU4OTYsMjguNDgwODU2OCAxNjcuNjA5Mjc1LDI4LjQyNzk2MjIgMTY3LjM0NjE2NCwyNy45ODIzMDQgYyAxNjcuMTI0MTk0LDI3LjYwNjMyOTMgMTY2LjYyMjc4LDI3LjE0NTk3NTkgMTY2LjI4NzY3MywyNi44NTU0MDg4IGMgMTY1Ljk0NDYzNSwyNi41NTc5NjQzIDE2NS43MTEwNTQsMjYuMjQ4NDk0NyAxNjUuNDU0MzMxLDI2LjAwNjQ1MzEgYyAxNjUuMTQ0NjI4LDI1LjcxNDQ2MTMgMTY1LjA0NTA2MSwyNS40ODQ1OTQgMTY0LjY3ODc2NiwyNS4yNzExMDM0IGMgMTY0LjM0MjI3LDI1LjA3NDk4MDUgMTYzLjkyMjM5OCwyNS4wNDY2NjE3IDE2My41NDI3OTgsMjQuODc3Mzk1NSBjIDE2My4xOTQ0NzgsMjQuNzIyMDc3MyAxNjMuMDgzOTcyLDI0LjM1NTIwNjIgMTYyLjUwNzUzNiwyNC4yMjgxODc3IGMgMTYxLjkzMTEsMjQuMTAxMTY5MyAxNjEuNDg2NTA1LDI0LjIyNzg5MzYgMTYxLjAyMzQ0MiwyNC4wNzI5NDE4IGMgMTYwLjY4MzkyNiwyMy45NTkzMzEzIDE2MC42MTM3MjIsMjMuNTYwMDE2OSAxNjAuMjQ3MjE0LDIzLjQzNzg0NzggYyAxNTkuNzg4MDMyLDIzLjI4NDc4NyAxNTkuMTQyODc5LDIzLjMwNzkwMzMgMTU4LjQ3MDEyMiwyMy4zNDc5MTE1IGMgMTU3Ljk0NDYxMiwyMy4zNzkxNjMgMTU3LjQwMjI2LDIzLjQyMDcyMTMgMTU2LjkxODU0NSwyMy4zOTY2NjgzIGMgMTU2LjU0MDUwNywyMy4zNzc4NzAxIDE1Ni4xOTgyODQsMjMuMzE4OTk3IDE1NS45Mjc5MDgsMjMuMTgzODA5MSBjIDE1NS43MDY2NTYsMjMuMDczMTgyOSAxNTUuNDkxMDYxLDIyLjYyNTQzMDggMTU0LjkxMjQxOCwyMi41NDg3MTMzIGMgMTU0LjMzMzc3NiwyMi40NzE5OTU3IDE1My44MTUyNjksMjIuNTQwODg5OSAxNTMuNDc3NzQzLDIyLjM5NTM3MzQgYyAxNTMuMTYyODU3LDIyLjI1OTYxNzcgMTUzLjA1OTc2MywyMS44NjY5NjEyIDE1Mi45MDc2OCwyMS44Mjg5NDAzIGMgMTUyLjYxOTg4NSwyMS43NTY5OTE1IDE1Mi4zNDQyOTYsMjEuNzA4OTgzNyAxNTIuMDc5NzExLDIxLjY3Nzc2NzEgYyAxNTEuNzIzNDg5LDIxLjYzNTczODYgMTUxLjM4NzIxMywyMS42MjQxNDY3IDE1MS4wNjc5NTEsMjEuNjI1NTQyMiBjIDE1MC41Mzk3MTIsMjEuNjI3ODUxMiAxNTAuMDU4MDUsMjEuNjY1NzE0MyAxNDkuNjA5Njc1LDIxLjY2MDA5NjEgYyAxNDkuMTY0MjU2LDIxLjY1NDUxNDkgMTQ4Ljc1MTY4OCwyMS42MDYwMjQgMTQ4LjM1ODk0MiwyMS40MzcxNDEzIGMgMTQ3Ljk3NTk4NywyMS4yNzI0NjkgMTQ3Ljg5ODgyNiwyMS4wOTUwNTI3IDE0Ny4yMDg2NSwyMC40MTc2MTgzIGMgMTQ2Ljk5NjQ3NiwyMC4yMDkzNjE0IDE0Ni43OTkxNzUsMjAuMDI2NDQ2NiAxNDYuNjIwNDgyLDE5Ljg2MjMxMzMgYyAxNDYuMjE3OTA5LDE5LjQ5MjU0MTMgMTQ1LjkwOTc4NCwxOS4yMTgwOTQ0IDE0NS43Mzg4MjYsMTguOTYzOTU2OCBjIDE0NS40OTE5ODQsMTguNTk3MDEzNCAxNDUuNTgzNTgxLDE4LjE0NTM5MDIgMTQ1LjM2ODE4MiwxNy43OTI1NTk3IGMgMTQ1LjE1Mjc4MywxNy40Mzk3MjkyIDE0NS4xNjU5MjYsMTcuNzk4NDk1MyAxNDQuODk2MDk0LDE2Ljk4ODk5OTUgYyAxNDQuODEwNTkzLDE2LjczMjQ5NTIgMTQ1LjI4NzU5OSwxNi42MDE0NDM0IDE0NS4zNjgxODIsMTYuMjk4ODIzNyBjIDE0NS42MzkzNzIsMTUuMjgwNDA2NyAxNDUuNjQzOTU3LDE1LjE4NTA5NjYgMTQ2LjE2NTk4MiwxNC42NTk0MjQ2IGMgMTQ2LjM2MDEyNywxNC40NjM5MjI4IDE0Ni40OTA2MjksMTMuNDk3Njc2NSAxNDYuNzI2MDksMTMuMjYyMjE2NCBjIDE0Ni45NTg3MzksMTMuMDI5NTY3NCAxNDcuMjIyNDM5LDEyLjgyMzA5MjEgMTQ3LjQ1OTI2OSwxMi41OTA1NTQzIGMgMTQ3LjY5MTkwOCwxMi4zNjIxMzE0IDE0Ny44OTg2MjEsMTIuMTA4NTYwMSAxNDguMDI0NTA2LDExLjc4MDMyODMgYyAxNDguMjc4NTQ1LDExLjExNzk1MTUgMTQ3Ljk4MjE2NywxMS4wNzkxMDUzIDE0OC4zNTg5NDIsMTAuOTA1MzA4OSBjIDE0OC43MzU3MTcsMTAuNzMxNTEyNSAxNDguNzcxMzA5LDExLjE2MDk0MiAxNDkuMjIxNjYxLDExLjM5OTI3MTYgYyAxNDkuNTA0NTY3LDExLjU0ODk4NzMgMTQ5LjgyMjMxNCwxMS41MjAzNjg3IDE1MC4xOTc5NDEsMTEuNTU0NTE2NiBjIDE1MC41NDY0MSwxMS41ODYxOTU2IDE1MS4xMjkyOTcsMTEuNTc3ODc4NyAxNTEuNzE2Mzk3LDExLjU2NzYyODEgYyAxNTIuNjA2MDI0LDExLjU1MjA5NTYgMTUzLjUwNTMyMiwxMS41MzIxMjM0IDE1My42MTMzNCwxMS42NDAxNDA5IGMgMTUzLjc2NjUwOCwxMS43OTMzMDg1IDE1My45NzE1NzcsMTIuMDg2NDY1NCAxNTQuMjA2MDk2LDEyLjIwMzcyNDkgYyAxNTQuNjk4MzE1LDEyLjQ0OTgzNDQgMTU1LjMwNjczLDEyLjQyNjY2MDggMTU1LjkyMDA3MywxMi4zODU1OTMzIGMgMTU2LjYxNDc1NSwxMi4zMzkwNzk3IDE1Ny4zMTU3NTksMTIuMjY5NjExNSAxNTcuODYxNDIsMTIuNTQyNDQyMSBjIDE1Ny45ODQ3MjgsMTIuNjA0MDk2MiAxNTguMjMzMzczLDEyLjk1ODM0OTMgMTU4LjM0MTQ2OCwxMy4wMTIzOTY5IGMgMTU5LjAyNDI4NCwxMy4zNTM4MDQ4IDE1OS4yNzgyMDQsMTMuMjYxNTc4OCAxNTkuNjM5Njg1LDEzLjQ1ODg3NTggYyAxNTkuOTQyMzA4LDEzLjYyNDA0NzMgMTYwLjExMzU1OSwxMy44OTY1MzY0IDE2MC40ODY0NzgsMTQuMjcxMDkyNyBjIDE2MC42ODY4MDEsMTQuNDcyMjk0NiAxNjAuODg2NjUyLDE0LjY4NTk4MSAxNjEuMTg0NywxNC45ODQwMjg1IGMgMTYxLjM2NTc4MywxNS4xNjUxMTEzIDE2MS42OTU3NjUsMTUuNTM3MDc2OCAxNjIuMDU4OTU3LDE1LjkzMjYyODcgYyAxNjIuMjY3MzksMTYuMTU5NjMzNiAxNjIuNDg2NzYxLDE2LjM5NDQwNjcgMTYyLjY5NTIwMiwxNi42MDUzMjY3IGMgMTYyLjk4MTUzLDE2Ljg5NTA2MDggMTYzLjI0NzIzNCwxNy4xMzk3ODQ4IDE2My40MzU2MzMsMTcuMjU3NTMzOCBjIDE2NC4wMDAxNjIsMTcuNjEwMzY0MyAxNjQuMDIxMjgyLDE3LjMyNTIwNSAxNjQuMzM5Mzg3LDE3LjQxMTUwMjYgYyAxNjQuODEwMjM4LDE3LjUzOTIzNzggMTY1LjI2NDA0NywxNy40NDE0NjIzIDE2NS42NTcxMjMsMTcuMTk4NDgwNSBjIDE2NS45OTU2NywxNi45ODkyMDYzIDE2Ni4yODkxNjYsMTYuNjcyMjE4OSAxNjYuNTA5Njk3LDE2LjI5ODgyMzcgYyAxNjYuNzAxOTU4LDE1Ljk3MzI5NSAxNjYuNzQxMDMsMTUuNDg4MDU1MSAxNjYuNzM0MTA0LDE0LjkxMzk5MjcgYyAxNjYuNzI4MDg0LDE0LjQxNDk5MzggMTY2LjY4NzMxLDEzLjg0ODg4MjIgMTY2LjY4MjE4MiwxMy4yNjIyMTY0IGMgMTY2LjY3OTAzNiwxMi45MDIzODQ4IDE2Ni42OTUxOTcsMTIuNTAwNjg5OCAxNjYuNjg3Njk1LDEyLjEyOTQzMjUgYyAxNjYuNjc3NDMyLDExLjYyMTU0MDUgMTY2LjYyMjg4MiwxMS4xNzA2MTMgMTY2LjQxNDAzMSwxMC45NjE3NjE3IGMgMTY2LjEyODAyOCwxMC42NzU3NTkyIDE2Ni4xNTkwNzYsMTAuNzkwNTcwOSAxNjUuOTYyNDA3LDEwLjM5NzIzMjkgYyAxNjUuODU2MjEyLDEwLjE4NDg0NDQgMTY1LjkyNDg3Nyw5LjU5MjQyMjA2IDE2NS40NTQzMzEsOS4xNjkzODI5NyBjIDE2NC45ODM3ODUsOC43NDYzNDM4OCAxNjQuOTc0NDgxLDguMzYxODA3OTcgMTY0LjQ1MjI5Miw4LjI2NjEzNjg1IGMgMTY0LjA2NjMwNyw4LjE5NTQyMDAyIDE2My42ODUyMjMsOC4xNDkwOTg1MiAxNjMuMzA2Mzk3LDguMTIwMTEyMDkgYyAxNjIuNzkxOTczLDguMDgwNzUwMjMgMTYyLjI4MTcxMSw4LjA3MzM1NDM1IDE2MS43Njg5OTQsOC4wODAyNDUyOCBjIDE2MS4yOTYwMjYsOC4wODY2MDE5NiAxNjAuODIwOTY5LDguMTA1MTE2MDcgMTYwLjMzODYyNiw4LjEyMTkwOTc3IGMgMTU5Ljg1MTQxNyw4LjEzODg3MjkxIDE1OS4zNTY3NzUsOC4xNTQwODA3NiAxNTguODQ5MzQ1LDguMTUzMjMxMTYgYyAxNTguMzA0NzEsOC4xNTIzMTkyNyAxNTcuNTIzNSw4LjE5NzkyNDc1IDE1Ny4wMjg3NCw4LjA3NDIzNDg4IGMgMTU2LjQ3ODY3Myw3LjkzNjcxNzk0IDE1Ni40MTA3NTcsNy41Mzc5NTM1MiAxNTUuOTI3OTA4LDcuMzc3MDA0MSBjIDE1NS40NDg2MjQsNy4yMTcyNDI1NSAxNTQuNjU5OTA5LDcuMzc3MDA0MSAxNTMuNzA3OCw3LjI3ODg3MzEzIGMgMTUzLjI0MTE2OSw3LjIzMDc3ODg1IDE1My4xNzc0MjQsNi45MTM2NTE5NCAxNTIuOTA3NjgsNi43Mjg0NTc2IGMgMTUyLjYzNzY4OCw2LjU0MzA5MzE0IDE1Mi4yMzA1NSw2LjQxNzk0MDg2IDE1MS41NjY5MjQsNi40NzQ0MTk2NiBjIDE1MS4xMjA3OCw2LjUxMjM4OTMzIDE1MC44MzIwNyw2LjU3MzIxMjE4IDE1MC4zNzE5NTMsNy4xMjM2Mjc3MSBjIDE0OS44MDI3NzEsNy4zNjM1NTI0MSAxNDkuMjc5NjkyLDcuMzcxMDU0MyAxNDkuMDEyNDMxLDcuNTA0Njg0NjEgYyAxNDguODc5MDg5LDcuOTQ5MTYxNTggMTQ4LjQzNDc0Niw4LjA5MzM0ODE2IDE0Ny45NTMwMTgsOC4xNTA3MDE5MyBjIDE0Ny4zNjQwMjMsOC4yMjA4MjY2MyAxNDYuNzE5MTM5LDguMTYxMTQzMjIgMTQ2LjUxODQ3NSw4LjM2MTgwNzk2IGMgMTQ2LjM5NDgyMSw4LjQ4NTQ2MTY5IDE0Ni4zMTQ5NTYsOC42NTM4MjIyNiAxNDYuMTY1OTgyLDguODAyNzk2NTkgYyAxNDUuOTQ4Njc4LDkuMDIwMTAwODcgMTQ1LjQ0MDE5MSw4Ljg5MDc5MTU4IDE0NS4xMzAzMDgsOS4wNDU3MzI3NiBjIDE0NC44ODA0Myw5LjE3MDY3MTk0IDE0NC42NjAyNiw5LjQ0ODYzNzQ2IDE0NC40NTYxMDcsOS42NTI3OTA3MiBjIDE0NC4zNDY4NDMsOS43NjIwNTQxMiAxNDMuNzk4MTcyLDkuNjcwODQ0OTIgMTQzLjMwNDMzNSwxMC4wMjk3MzE2IGMgMTQyLjgxMDQ5OCwxMC4zODg2MTg0IDE0Mi4xNTc4NjUsMTEuMDQzMjUxMiAxNDIuMDg5ODUsMTEuMTM1NTU4MiBjIDE0MS45NTc2MjgsMTEuMzE1MDAxMSAxNDEuNjQzNDI0LDExLjUxMDg3NDMgMTQxLjU0MzksMTEuNzA5OTIyNCBjIDE0MS4xMDI3NTIsMTIuNTkyMjE3NiAxNDEuMzU3ODIyLDEyLjQ3MDM4MDkgMTQxLjIzNjQ5MywxMi42OTgwMTQ0IGMgMTQxLjE3NjQwNSwxMi44MTA3NTA1IDE0MS4wNjEzNTgsMTMuMDQ3NTM5NiAxNDAuNjUwMDIyLDEzLjQ1ODg3NTkgYyAxNDAuNjQyNDU0LDEzLjQ2NjQ0MzkgMTQwLjQ3NTMyNiwxMy44NjcxNDQ1IDE0MC40Nzk0NDEsMTQuMjcxMDkyNCBjIDE0MC40ODMyMzIsMTQuNjQzMzM5NSAxNDAuNDk2MjQ0LDE1LjAwNTU4MzIgMTQwLjQwOTQyNCwxNS4xNzkyMjM1IGMgMTQwLjExNTQwMSwxNS43NjcyNjg2IDE0MC4wNzkyOTIsMTUuNTI4ODk5NyAxMzkuODQ4MDA3LDE1Ljk5MzUzNDQgYyAxMzkuNzU3OTQxLDE2LjE3NDQ3MDEgMTM5LjU4MDM2LDE2LjY4ODkyNzggMTM5LjY5OTMyNCwxNi45ODg5OTkyIGMgMTM5Ljg4NTg1NCwxNy40NTk0OTc5IDE0MC4zNzQ0NzYsMTcuOTAyODg4NiAxNDAuNDA5NDI0LDE4LjA3NzYyNzQgYyAxNDAuNDc2NjgsMTguNDEzOTA4IDE0MC40NzkzMzQsMTguNjYyODMzIDE0MC40Nzk0NDEsMTkuMTE3MDM5OCBjIDE0MC40Nzk1NjksMTkuNjY2NTk3OSAxNDAuNzM5MDE5LDE5Ljk2MjEzODMgMTQxLjA1MzM3LDIwLjIyNjYyNDUgYyAxNDEuMzQ0Njg2LDIwLjQ3MTcyODggMTQxLjY4MzE1MSwyMC42OTAxNjM1IDE0MS45MDYwNzMsMjEuMDU5Mzc5MiBjIDE0Mi4zNjk1NDcsMjEuODI3MDA2NyAxNDIuMjA1ODg0LDIyLjMzMTE5ODYgMTQyLjM2OTU0NywyMi4zOTUzNzM0IGMgMTQyLjc2MDYwMywyMi41NDg3MTMzIDE0My4zMjg4NjYsMjIuNDk2NjU1NiAxNDMuNjAwNjQ5LDIyLjc2ODQzOCBjIDE0My44ODM4MjEsMjMuMDUxNjEgMTQzLjgzMjM4NSwyMy45OTk1NDAzIDE0NC4wOTMwODksMjQuMDcyOTQxOCBjIDE0NC43ODY3MywyNC4yNjgyMzY5IDE0NC45MTg0NSwyNC4zMzI2NjA0IDE0NS4yNjY0MjYsMjQuNDQ4NTI5NSBjIDE0NS42MTQ0MDMsMjQuNTY0Mzk4NyAxNDUuNjg4MDI1LDI0LjkyNjQ4NzkgMTQ2LjE2NTk4MiwyNS4zMzIwMjYyIGMgMTQ2LjY0MzkzOCwyNS43Mzc1NjQ2IDE0Ni44ODg1ODQsMjYuMDA2NDUzMSAxNDcuMjA4NjUsMjYuMzc0ODQxMiBjIDE0Ny40NzM1NzYsMjYuNjc5NzYzNCAxNDcuODI4OTYsMjYuNzAyMjYyNiAxNDguMTc4NTEzLDI2Ljc5MjY1MDcgYyAxNDguNDUzMjgzLDI2Ljg2MzcwMSAxNDguNzI0NDQ5LDI2Ljk3NjY5ODkgMTQ4Ljk0NTI0NywyNy4zMDE3ODg4IGMgMTQ5LjEwNTI4MSwyNy41Mzc0MTM3IDE0OS41OTIyNCwyNy41ODk3NjIyIDE1MC4xNzQ3MDgsMjcuNTkxMDQxMSBjIDE1MC42MjU2ODksMjcuNTkyMDMxMyAxNTEuMTMzOTI3LDI3LjU2MjQwNjIgMTUxLjU5MjAwNywyNy41NjM1MzAyIGMgMTUyLjA1MTY4LDI3LjU2NDY1ODIgMTUyLjQ2MDg0NywyNy41OTY3NDk0IDE1Mi43MTA5NywyNy43MjE4MTA4IGMgMTUzLjExOTc3NiwyNy45MjYyMTM4IDE1Mi45NDc4NiwyNy45MDQ1MjQ1IDE1My4xOTAxNjksMjguMTQ2ODMzNCBjIDE1My41MDM2ODIsMjguNDYwMzQ2IDE1NC4wNzM2NTksMjguNDU2OTgzOSAxNTQuNjU0MDYyLDI4LjQzNTczNTMgYyAxNTUuMTg1NDU2LDI4LjQxNjI4MSAxNTUuNzI1NTksMjguMzgxODMzMyAxNTYuMDg1NjM2LDI4LjU2MTg1NjYgYyAxNTYuNDM3OTY2LDI4LjczODAyMTEgMTU2LjM2OTgwNSwyOC45Mjg0MzYzIDE1Ni45NDAxNjYsMjkuMjEzNjE2OCBjIDE1Ny4wNzQzMjQsMjkuMjgwNjk2IDE1Ny40OTcyNjQsMjkuMjY3Njk3IDE1Ny45MDg4NDgsMjkuMjkyNTE4MyBjIDE1OC4xODY4MDgsMjkuMzA5MjgxMiAxNTguNDU5NTg4LDI5LjM0MzI5MzUgMTU4LjYzNDc0MSwyOS40MzA4Njk2IGMgMTU5LjA2OTI0OCwyOS42NDgxMjI0IDE1OS4wMzI1NSwyOS44NTk4ODk3IDE1OS4xOTU3NTYsMjkuOTE0MjkxOSBjIDE1OS41Nzk3OTcsMzAuMDQyMzA1NCAxNjAuMzg2NjE2LDMwLjE2Njk2NTQgMTYwLjYwOTU1NCwzMC4zODk5MDMyIGMgMTYwLjk5ODk2LDMwLjc3OTMwOTUgMTYxLjM1OTE4OCwzMS4yNTA4NjI4IDE2MS42NDczMTksMzEuNTE2NTAxMyBjIDE2MS45ODc1ODIsMzEuODMwMjAzMSAxNjIuMjEyMzgsMzIuMDYxNzU5IDE2Mi4zMzI1MzcsMzIuNTEzNTAzIGMgMTYyLjQwMzcwMywzMi43ODEwNjIzIDE2Mi40MzgxNjIsMzMuMTI1ODYyOSAxNjIuNDM4MTYyLDMzLjYxMDcyMDcgYyAxNjIuNDM4MTYyLDM0LjAxODY2MTkgMTYyLjQ1MjgzNywzNC43Mzc0MDg1IDE2Mi4yMjY2NjEsMzQuOTYzNTg0NSBjIDE2Mi4wOTk0MTUsMzUuMDkwODI5NiAxNjEuMzQzMDQsMzUuMTc0MDg0OCAxNjEuMTg0Njk5LDM1LjI1MzI1NTUgYyAxNjAuODAzMTUzLDM1LjQ0NDAyODcgMTYwLjQzNjYzNiwzNS44OTkwODE1IDE2MC4xMjY1NDcsMzYuMjA5MTcwNCBjIDE1OS44NDE3OTQsMzYuNDkzOTIzMyAxNTkuMTc5MTc3LDM3LjI0NDIyNTEgMTU4Ljg0OTM0NCwzNy41MjcxNzI3IGMgMTU4LjU2MjMyMywzNy43NzMzOTQ3IDE1Ny45NDU4MywzNy42ODg3ODMxIDE1Ny42NDk4NTksMzcuODQ1ODEyOSBjIDE1Ny4yODg1NTksMzguMDM3NTAzNSAxNTcuMzA0MzM0LDM4LjIyOTQyMzIgMTU3LjAyODczOSwzOC4zNjcyMjAzIGMgMTU2Ljc3OTYxNSwzOC40OTE3ODI3IDE1Ni4zNDkyMzIsMzguNTMyNTQ0NSAxNTUuODY5OTE1LDM4LjU0NDQ2ODEgYyAxNTUuNDgyNDQxLDM4LjU1NDEwNyAxNTUuMDYyOTg4LDM4LjU0NDkwMDUgMTU0LjY4MTQ2LDM4LjU0NTg4MzggYyAxNTQuMTMxODksMzguNTQ3MzAwMiAxNTMuNjYxMDEsMzguNTY5ODU5MyAxNTMuNDc3NzQyLDM4LjcwMDM0MDUgYyAxNTMuMDg5MTcyLDM4Ljk3Njk5MDQgMTUzLjAxOTIwMiwzOS4yNDYyOTIgMTUyLjQ0ODcyNSwzOS4zMzUwNDI0IGMgMTUyLjE4MTUzOSwzOS4zNzY2MDkxIDE1MS44MDQ1NjQsMzkuMzc4NTcwOCAxNTEuMjMzNjQ4LDM5LjMyMzEzMzQgYyAxNTAuNDE4NzQxLDM5LjI0NDAwMzkgMTUwLjU5NjEwNiwzOC44NTk2NjA1IDE1MC4xOTc5NDEsMzguNzAwMzQwNSBjIDE1MC4wMjkwMDUsMzguNjMyNzQzIDE0OS42ODAzODQsMzguNTk1Mzc0MSAxNDkuMzIxMDMzLDM4LjU2MDgxOTggYyAxNDguODMzNDMzLDM4LjUxMzkzMzQgMTQ4LjMyNjA3OSwzOC40NzIyMjkyIDE0OC4yMjEwNywzOC4zNjcyMjAzIGMgMTQ3Ljk2MjgxNCwzOC4xMDg5NjQ2IDE0OC4wNDU4OTgsMzguMDc2NzI5OSAxNDcuNzg2NTYzLDM3Ljk1NDgzNzUgYyAxNDcuNDc1NTI3LDM3LjgwODY0NDEgMTQ2Ljk3Nzk0NSwzNy43MTI1NTE3IDE0Ni42MDY3MjksMzcuNDY3NzA3NCBjIDE0Ni4yMTg3OTYsMzcuMjExODM3MiAxNDUuODQ3MTQsMzYuODEwNjQ4MiAxNDUuNDY1MTY3LDM2LjQwMjk2MzIgYyAxNDUuMTc2ODcxLDM2LjA5NTI2MTMgMTQ0Ljg4MjY5OSwzNS43ODM4NTg5IDE0NC41NzEyMTYsMzUuNTI4NDQyOSBjIDE0NC4zMzAzODcsMzUuMzMwOTYyNyAxNDMuNzMxMTcsMzUuMTY5MDkyOSAxNDMuMzA0MzM1LDM1LjAzNjAwMjcgYyAxNDIuODc3NSwzNC45MDI5MTI1IDE0My4wODA5LDM0LjA1Nzc3ODQgMTQyLjg3MzkzNiwzMy44MTU0MjcxIGMgMTQyLjYwMDY2NywzMy40OTU0MzM4IDE0Mi4xNDgzNTksMzMuMzc5MDIzMyAxNDEuNDQ2OTE5LDMzLjQxMzg4MTUgYyAxNDEuMDQ3NzQsMzMuNDMzNzE4OCAxNDAuNTY3ODc4LDMzLjUwMjU0NTggMTM5Ljk5NDQxNSwzMy42MTA3MjA3IGMgMTM5LjYyNDg3LDMzLjY4MDQyOTYgMTM5LjYxNzg0MywzNC4wOTQ1NzE0IDEzOC45NTE1OTksMzQuNTg4NDY5MSBsIDEzOC44MTE1MTcsMzUuMTA2MTAxMSB6IFwiO1xuLy8gdmFyIHBhdGggPSBwb2ludDEgKyBwb2ludDIgKyBwb2ludDMgKyBwb2ludDQgKyBwb2ludDUgKyBwb2ludDYgKyBwb2ludDcgKyBwb2ludDggKyBwb2ludDk7XG52YXIgcGF0aCA9IGxldHRyZXNBTlQgKyBsZXR0cmVTO1xuXG5tb2R1bGUuZXhwb3J0cyA9IHBhdGg7IiwiJ3VzZSBzdHJpY3QnO1xuXG52YXIgc3FydCA9IE1hdGguc3FydDtcbnZhciBwb3cgPSBNYXRoLnBvdztcblxuZnVuY3Rpb24gc2lnbih4KSB7XG5cdHJldHVybiB4ID8geCA8IDAgPyAtMSA6IDEgOiAwO1xufVxuXG5mdW5jdGlvbiByYW5nZShzdGFydCwgY291bnQpIHtcbiAgICByZXR1cm4gQXJyYXkuYXBwbHkoMCwgQXJyYXkoY291bnQpKS5tYXAoZnVuY3Rpb24gKGVsZW1lbnQsIGluZGV4KSB7XG4gICAgXHRyZXR1cm4gaW5kZXggKyBzdGFydFxuICAgIH0pO1xufVxuXG5mdW5jdGlvbiBkaXN0YW5jZShhLCBiKXtcblx0cmV0dXJuIHNxcnQocG93KGEueCAtIGIueCwgMikgKyBwb3coYS55IC0gYi55LCAyKSk7XG59XG5cbmZ1bmN0aW9uIG5vcm0odil7XG5cdHJldHVybiBzcXJ0KHBvdyh2LngsIDIpICsgcG93KHYueSwgMikpO1xufVxuXG5tb2R1bGUuZXhwb3J0cyA9IHtcblx0c2lnbjogc2lnbixcblx0cmFuZ2U6IHJhbmdlLFxuXHRkaXN0YW5jZTogZGlzdGFuY2UsXG5cdG5vcm06IG5vcm1cbn0iLCIndXNlIHN0cmljdCdcblxuZnVuY3Rpb24gVmVjdG9yKHgsIHkpIHtcbiAgICB0aGlzLnggPSB4OyAgICAgICAgICAgICAgICBcbiAgICB0aGlzLnkgPSB5O1xufVxuXG5WZWN0b3IucHJvdG90eXBlLm5vcm0gPSBmdW5jdGlvbigpe1xuXHRyZXR1cm4gTWF0aC5zcXJ0KHRoaXMueCAqIHRoaXMueCArIHRoaXMueSAqIHRoaXMueSk7XG59XG5cblZlY3Rvci5wcm90b3R5cGUubm9ybWFsaXplID0gZnVuY3Rpb24oKXtcblx0dmFyIG5vcm0gPSB0aGlzLm5vcm0oKTtcblx0dGhpcy54ID0gdGhpcy54IC8gbm9ybTtcblx0dGhpcy55ID0gdGhpcy55IC8gbm9ybTtcbn1cblxuXG5cbm1vZHVsZS5leHBvcnRzID0gVmVjdG9yOyJdfQ==
