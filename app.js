(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js":[function(require,module,exports){
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
},{"robust-orientation":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js","simplicial-complex":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js"}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js":[function(require,module,exports){
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
},{}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js":[function(require,module,exports){
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
},{"two-product":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js","two-sum":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/node_modules/two-sum/two-sum.js"}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js":[function(require,module,exports){
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
},{}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js":[function(require,module,exports){
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
},{}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js":[function(require,module,exports){
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
},{}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/orientation.js":[function(require,module,exports){
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
},{"robust-scale":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-scale/robust-scale.js","robust-subtract":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-subtract/robust-diff.js","robust-sum":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/robust-sum/robust-sum.js","two-product":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/robust-orientation/node_modules/two-product/two-product.js"}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js":[function(require,module,exports){
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


},{}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js":[function(require,module,exports){
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
},{}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/topology.js":[function(require,module,exports){
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

},{"bit-twiddle":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/bit-twiddle/twiddle.js","union-find":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/node_modules/simplicial-complex/node_modules/union-find/index.js"}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js":[function(require,module,exports){
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

},{}],"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/triangulate.js":[function(require,module,exports){
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
},{"incremental-convex-hull":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/incremental-convex-hull/ich.js","uniq":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/node_modules/uniq/uniq.js"}],"/Users/vallette/Desktop/colony/node_modules/game-shell/lib/hrtime-polyfill.js":[function(require,module,exports){
if(typeof window.performance === "object") {
  if(window.performance.now) {
    module.exports = function() { return window.performance.now() }
  } else if(window.performance.webkitNow) {
    module.exports = function() { return window.performance.webkitNow() }
  }
} else if(Date.now) {
  module.exports = Date.now
} else {
  module.exports = function() { return (new Date()).getTime() }
}

},{}],"/Users/vallette/Desktop/colony/node_modules/game-shell/lib/mousewheel-polyfill.js":[function(require,module,exports){
//Adapted from here: https://developer.mozilla.org/en-US/docs/Web/Reference/Events/wheel?redirectlocale=en-US&redirectslug=DOM%2FMozilla_event_reference%2Fwheel

var prefix = "", _addEventListener, onwheel, support;

// detect event model
if ( window.addEventListener ) {
  _addEventListener = "addEventListener";
} else {
  _addEventListener = "attachEvent";
  prefix = "on";
}

// detect available wheel event
support = "onwheel" in document.createElement("div") ? "wheel" : // Modern browsers support "wheel"
          document.onmousewheel !== undefined ? "mousewheel" : // Webkit and IE support at least "mousewheel"
          "DOMMouseScroll"; // let's assume that remaining browsers are older Firefox

function _addWheelListener( elem, eventName, callback, useCapture ) {
  elem[ _addEventListener ]( prefix + eventName, support == "wheel" ? callback : function( originalEvent ) {
    !originalEvent && ( originalEvent = window.event );

    // create a normalized event object
    var event = {
      // keep a ref to the original event object
      originalEvent: originalEvent,
      target: originalEvent.target || originalEvent.srcElement,
      type: "wheel",
      deltaMode: originalEvent.type == "MozMousePixelScroll" ? 0 : 1,
      deltaX: 0,
      delatZ: 0,
      preventDefault: function() {
        originalEvent.preventDefault ?
          originalEvent.preventDefault() :
          originalEvent.returnValue = false;
      }
    };
    
    // calculate deltaY (and deltaX) according to the event
    if ( support == "mousewheel" ) {
      event.deltaY = - 1/40 * originalEvent.wheelDelta;
      // Webkit also support wheelDeltaX
      originalEvent.wheelDeltaX && ( event.deltaX = - 1/40 * originalEvent.wheelDeltaX );
    } else {
      event.deltaY = originalEvent.detail;
    }

    // it's time to fire the callback
    return callback( event );
  }, useCapture || false );
}

module.exports = function( elem, callback, useCapture ) {
  _addWheelListener( elem, support, callback, useCapture );

  // handle MozMousePixelScroll in older Firefox
  if( support == "DOMMouseScroll" ) {
    _addWheelListener( elem, "MozMousePixelScroll", callback, useCapture );
  }
};
},{}],"/Users/vallette/Desktop/colony/node_modules/game-shell/lib/raf-polyfill.js":[function(require,module,exports){
// http://paulirish.com/2011/requestanimationframe-for-smart-animating/
// http://my.opera.com/emoller/blog/2011/12/20/requestanimationframe-for-smart-er-animating
 
// requestAnimationFrame polyfill by Erik Möller. fixes from Paul Irish and Tino Zijdel
 
// MIT license
var lastTime = 0;
var vendors = ['ms', 'moz', 'webkit', 'o'];
for(var x = 0; x < vendors.length && !window.requestAnimationFrame; ++x) {
    window.requestAnimationFrame = window[vendors[x]+'RequestAnimationFrame'];
    window.cancelAnimationFrame = window[vendors[x]+'CancelAnimationFrame'] 
                               || window[vendors[x]+'CancelRequestAnimationFrame'];
}

if (!window.requestAnimationFrame)
    window.requestAnimationFrame = function(callback, element) {
        var currTime = new Date().getTime();
        var timeToCall = Math.max(0, 16 - (currTime - lastTime));
        var id = window.setTimeout(function() { callback(currTime + timeToCall); }, 
          timeToCall);
        lastTime = currTime + timeToCall;
        return id;
    };

if (!window.cancelAnimationFrame)
    window.cancelAnimationFrame = function(id) {
        clearTimeout(id);
    };

},{}],"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/binary-search-bounds/search-bounds.js":[function(require,module,exports){
"use strict"

function compileSearch(funcName, predicate, reversed, extraArgs, useNdarray, earlyOut) {
  var code = [
    "function ", funcName, "(a,l,h,", extraArgs.join(","),  "){",
earlyOut ? "" : "var i=", (reversed ? "l-1" : "h+1"),
";while(l<=h){\
var m=(l+h)>>>1,x=a", useNdarray ? ".get(m)" : "[m]"]
  if(earlyOut) {
    if(predicate.indexOf("c") < 0) {
      code.push(";if(x===y){return m}else if(x<=y){")
    } else {
      code.push(";var p=c(x,y);if(p===0){return m}else if(p<=0){")
    }
  } else {
    code.push(";if(", predicate, "){i=m;")
  }
  if(reversed) {
    code.push("l=m+1}else{h=m-1}")
  } else {
    code.push("h=m-1}else{l=m+1}")
  }
  code.push("}")
  if(earlyOut) {
    code.push("return -1};")
  } else {
    code.push("return i};")
  }
  return code.join("")
}

function compileBoundsSearch(predicate, reversed, suffix, earlyOut) {
  var result = new Function([
  compileSearch("A", "x" + predicate + "y", reversed, ["y"], false, earlyOut),
  compileSearch("B", "x" + predicate + "y", reversed, ["y"], true, earlyOut),
  compileSearch("P", "c(x,y)" + predicate + "0", reversed, ["y", "c"], false, earlyOut),
  compileSearch("Q", "c(x,y)" + predicate + "0", reversed, ["y", "c"], true, earlyOut),
"function dispatchBsearch", suffix, "(a,y,c,l,h){\
if(a.shape){\
if(typeof(c)==='function'){\
return Q(a,(l===undefined)?0:l|0,(h===undefined)?a.shape[0]-1:h|0,y,c)\
}else{\
return B(a,(c===undefined)?0:c|0,(l===undefined)?a.shape[0]-1:l|0,y)\
}}else{\
if(typeof(c)==='function'){\
return P(a,(l===undefined)?0:l|0,(h===undefined)?a.length-1:h|0,y,c)\
}else{\
return A(a,(c===undefined)?0:c|0,(l===undefined)?a.length-1:l|0,y)\
}}}\
return dispatchBsearch", suffix].join(""))
  return result()
}

module.exports = {
  ge: compileBoundsSearch(">=", false, "GE"),
  gt: compileBoundsSearch(">", false, "GT"),
  lt: compileBoundsSearch("<", true, "LT"),
  le: compileBoundsSearch("<=", true, "LE"),
  eq: compileBoundsSearch("-", true, "EQ", true)
}

},{}],"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/domready/ready.js":[function(require,module,exports){
/*!
  * domready (c) Dustin Diaz 2014 - License MIT
  */
!function (name, definition) {

  if (typeof module != 'undefined') module.exports = definition()
  else if (typeof define == 'function' && typeof define.amd == 'object') define(definition)
  else this[name] = definition()

}('domready', function () {

  var fns = [], listener
    , doc = document
    , domContentLoaded = 'DOMContentLoaded'
    , loaded = /^loaded|^c/.test(doc.readyState)

  if (!loaded)
  doc.addEventListener(domContentLoaded, listener = function () {
    doc.removeEventListener(domContentLoaded, listener)
    loaded = 1
    while (listener = fns.shift()) listener()
  })

  return function (fn) {
    loaded ? fn() : fns.push(fn)
  }

});

},{}],"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/invert-hash/invert.js":[function(require,module,exports){
"use strict"

function invert(hash) {
  var result = {}
  for(var i in hash) {
    if(hash.hasOwnProperty(i)) {
      result[hash[i]] = i
    }
  }
  return result
}

module.exports = invert
},{}],"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/iota-array/iota.js":[function(require,module,exports){
"use strict"

function iota(n) {
  var result = new Array(n)
  for(var i=0; i<n; ++i) {
    result[i] = i
  }
  return result
}

module.exports = iota
},{}],"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/uniq/uniq.js":[function(require,module,exports){
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

},{}],"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/vkey/index.js":[function(require,module,exports){
var ua = typeof window !== 'undefined' ? window.navigator.userAgent : ''
  , isOSX = /OS X/.test(ua)
  , isOpera = /Opera/.test(ua)
  , maybeFirefox = !/like Gecko/.test(ua) && !isOpera

var i, output = module.exports = {
  0:  isOSX ? '<menu>' : '<UNK>'
, 1:  '<mouse 1>'
, 2:  '<mouse 2>'
, 3:  '<break>'
, 4:  '<mouse 3>'
, 5:  '<mouse 4>'
, 6:  '<mouse 5>'
, 8:  '<backspace>'
, 9:  '<tab>'
, 12: '<clear>'
, 13: '<enter>'
, 16: '<shift>'
, 17: '<control>'
, 18: '<alt>'
, 19: '<pause>'
, 20: '<caps-lock>'
, 21: '<ime-hangul>'
, 23: '<ime-junja>'
, 24: '<ime-final>'
, 25: '<ime-kanji>'
, 27: '<escape>'
, 28: '<ime-convert>'
, 29: '<ime-nonconvert>'
, 30: '<ime-accept>'
, 31: '<ime-mode-change>'
, 27: '<escape>'
, 32: '<space>'
, 33: '<page-up>'
, 34: '<page-down>'
, 35: '<end>'
, 36: '<home>'
, 37: '<left>'
, 38: '<up>'
, 39: '<right>'
, 40: '<down>'
, 41: '<select>'
, 42: '<print>'
, 43: '<execute>'
, 44: '<snapshot>'
, 45: '<insert>'
, 46: '<delete>'
, 47: '<help>'
, 91: '<meta>'  // meta-left -- no one handles left and right properly, so we coerce into one.
, 92: '<meta>'  // meta-right
, 93: isOSX ? '<meta>' : '<menu>'      // chrome,opera,safari all report this for meta-right (osx mbp).
, 95: '<sleep>'
, 106: '<num-*>'
, 107: '<num-+>'
, 108: '<num-enter>'
, 109: '<num-->'
, 110: '<num-.>'
, 111: '<num-/>'
, 144: '<num-lock>'
, 145: '<scroll-lock>'
, 160: '<shift-left>'
, 161: '<shift-right>'
, 162: '<control-left>'
, 163: '<control-right>'
, 164: '<alt-left>'
, 165: '<alt-right>'
, 166: '<browser-back>'
, 167: '<browser-forward>'
, 168: '<browser-refresh>'
, 169: '<browser-stop>'
, 170: '<browser-search>'
, 171: '<browser-favorites>'
, 172: '<browser-home>'

  // ff/osx reports '<volume-mute>' for '-'
, 173: isOSX && maybeFirefox ? '-' : '<volume-mute>'
, 174: '<volume-down>'
, 175: '<volume-up>'
, 176: '<next-track>'
, 177: '<prev-track>'
, 178: '<stop>'
, 179: '<play-pause>'
, 180: '<launch-mail>'
, 181: '<launch-media-select>'
, 182: '<launch-app 1>'
, 183: '<launch-app 2>'
, 186: ';'
, 187: '='
, 188: ','
, 189: '-'
, 190: '.'
, 191: '/'
, 192: '`'
, 219: '['
, 220: '\\'
, 221: ']'
, 222: "'"
, 223: '<meta>'
, 224: '<meta>'       // firefox reports meta here.
, 226: '<alt-gr>'
, 229: '<ime-process>'
, 231: isOpera ? '`' : '<unicode>'
, 246: '<attention>'
, 247: '<crsel>'
, 248: '<exsel>'
, 249: '<erase-eof>'
, 250: '<play>'
, 251: '<zoom>'
, 252: '<no-name>'
, 253: '<pa-1>'
, 254: '<clear>'
}

for(i = 58; i < 65; ++i) {
  output[i] = String.fromCharCode(i)
}

// 0-9
for(i = 48; i < 58; ++i) {
  output[i] = (i - 48)+''
}

// A-Z
for(i = 65; i < 91; ++i) {
  output[i] = String.fromCharCode(i)
}

// num0-9
for(i = 96; i < 106; ++i) {
  output[i] = '<num-'+(i - 96)+'>'
}

// F1-F24
for(i = 112; i < 136; ++i) {
  output[i] = 'F'+(i-111)
}

},{}],"/Users/vallette/Desktop/colony/node_modules/game-shell/shell.js":[function(require,module,exports){
"use strict"

var EventEmitter = require("events").EventEmitter
  , util         = require("util")
  , domready     = require("domready")
  , vkey         = require("vkey")
  , invert       = require("invert-hash")
  , uniq         = require("uniq")
  , bsearch      = require("binary-search-bounds")
  , iota         = require("iota-array")
  , min          = Math.min

//Browser compatibility hacks
require("./lib/raf-polyfill.js")
var addMouseWheel = require("./lib/mousewheel-polyfill.js")
var hrtime = require("./lib/hrtime-polyfill.js")

//Remove angle braces and other useless crap
var filtered_vkey = (function() {
  var result = new Array(256)
    , i, j, k
  for(i=0; i<256; ++i) {
    result[i] = "UNK"
  }
  for(i in vkey) {
    k = vkey[i]
    if(k.charAt(0) === '<' && k.charAt(k.length-1) === '>') {
      k = k.substring(1, k.length-1)
    }
    k = k.replace(/\s/g, "-")
    result[parseInt(i)] = k
  }
  return result
})()

//Compute minimal common set of keyboard functions
var keyNames = uniq(Object.keys(invert(filtered_vkey)))

//Translates a virtual keycode to a normalized keycode
function virtualKeyCode(key) {
  return bsearch.eq(keyNames, key)
}

//Maps a physical keycode to a normalized keycode
function physicalKeyCode(key) {
  return virtualKeyCode(filtered_vkey[key])
}

//Game shell
function GameShell() {
  EventEmitter.call(this)
  this._curKeyState  = new Array(keyNames.length)
  this._pressCount   = new Array(keyNames.length)
  this._releaseCount = new Array(keyNames.length)
  
  this._tickInterval = null
  this._rafHandle = null
  this._tickRate = 0
  this._lastTick = hrtime()
  this._frameTime = 0.0
  this._paused = true
  this._width = 0
  this._height = 0
  
  this._wantFullscreen = false
  this._wantPointerLock = false
  this._fullscreenActive = false
  this._pointerLockActive = false
  
  this._render = render.bind(undefined, this)
  
  for(var i=0; i<keyNames.length; ++i) {
    this._curKeyState[i] = false
    this._pressCount[i] = this._releaseCount[i] = 0
  }
  
  //Public members
  this.element = null
  this.bindings = {}
  this.frameSkip = 100.0
  this.tickCount = 0
  this.frameCount = 0
  this.startTime = hrtime()
  this.tickTime = this._tickRate
  this.frameTime = 10.0
  this.stickyFullscreen = false
  this.stickyPointLock = false
  
  //Scroll stuff
  this.scroll = [0,0,0]
    
  //Mouse state
  this.mouseX = 0
  this.mouseY = 0
  this.prevMouseX = 0
  this.prevMouseY = 0
}

util.inherits(GameShell, EventEmitter)

var proto = GameShell.prototype

//Bind keynames
proto.keyNames = keyNames

//Binds a virtual keyboard event to a physical key
proto.bind = function(virtual_key) {
  //Look up previous key bindings
  var arr
  if(virtual_key in this.bindings) {
    arr = this.bindings[virtual_key]
  } else {
    arr = []
  }
  //Add keys to list
  var physical_key
  for(var i=1, n=arguments.length; i<n; ++i) {
    physical_key = arguments[i]
    if(virtualKeyCode(physical_key) >= 0) {
      arr.push(physical_key)
    } else if(physical_key in this.bindings) {
      var keybinds = this.bindings[physical_key]
      for(var j=0; j<keybinds.length; ++j) {
        arr.push(keybinds[j])
      }
    }
  }
  //Remove any duplicate keys
  arr = uniq(arr)
  if(arr.length > 0) {
    this.bindings[virtual_key] = arr
  }
  this.emit('bind', virtual_key, arr)
}

//Unbinds a virtual keyboard event
proto.unbind = function(virtual_key) {
  if(virtual_key in this.bindings) {
    delete this.bindings[virtual_key]
  }
  this.emit('unbind', virtual_key)
}

//Checks if a key is set in a given state
function lookupKey(state, bindings, key) {
  if(key in bindings) {
    var arr = bindings[key]
    for(var i=0, n=arr.length; i<n; ++i) {
      if(state[virtualKeyCode(arr[i])]) {
        return true
      }
    }
    return false
  }
  var kc = virtualKeyCode(key)
  if(kc >= 0) {
    return state[kc]
  }
  return false
}

//Checks if a key is set in a given state
function lookupCount(state, bindings, key) {
  if(key in bindings) {
    var arr = bindings[key], r = 0
    for(var i=0, n=arr.length; i<n; ++i) {
      r += state[virtualKeyCode(arr[i])]
    }
    return r
  }
  var kc = virtualKeyCode(key)
  if(kc >= 0) {
    return state[kc]
  }
  return 0
}

//Checks if a key (either physical or virtual) is currently held down
proto.down = function(key) {
  return lookupKey(this._curKeyState, this.bindings, key)
}

//Checks if a key was ever down
proto.wasDown = function(key) {
  return this.down(key) || !!this.press(key)
}

//Opposite of down
proto.up = function(key) {
  return !this.down(key)
}

//Checks if a key was released during previous frame
proto.wasUp = function(key) {
  return this.up(key) || !!this.release(key)
}

//Returns the number of times a key was pressed since last tick
proto.press = function(key) {
  return lookupCount(this._pressCount, this.bindings, key)
}

//Returns the number of times a key was released since last tick
proto.release = function(key) {
  return lookupCount(this._releaseCount, this.bindings, key)
}

//Pause/unpause the game loop
Object.defineProperty(proto, "paused", {
  get: function() {
    return this._paused
  },
  set: function(state) {
    var ns = !!state
    if(ns !== this._paused) {
      if(!this._paused) {
        this._paused = true
        this._frameTime = min(1.0, (hrtime() - this._lastTick) / this._tickRate)
        clearInterval(this._tickInterval)
        //cancelAnimationFrame(this._rafHandle)
      } else {
        this._paused = false
        this._lastTick = hrtime() - Math.floor(this._frameTime * this._tickRate)
        this._tickInterval = setInterval(tick, this._tickRate, this)
        this._rafHandle = requestAnimationFrame(this._render)
      }
    }
  }
})

//Fullscreen state toggle

function tryFullscreen(shell) {
  //Request full screen
  var elem = shell.element
  
  if(shell._wantFullscreen && !shell._fullscreenActive) {
    var fs = elem.requestFullscreen ||
             elem.requestFullScreen ||
             elem.webkitRequestFullscreen ||
             elem.webkitRequestFullScreen ||
             elem.mozRequestFullscreen ||
             elem.mozRequestFullScreen ||
             function() {}
    fs.call(elem)
  }
  if(shell._wantPointerLock && !shell._pointerLockActive) {
    var pl =  elem.requestPointerLock ||
              elem.webkitRequestPointerLock ||
              elem.mozRequestPointerLock ||
              elem.msRequestPointerLock ||
              elem.oRequestPointerLock ||
              function() {}
    pl.call(elem)
  }
}

var cancelFullscreen = document.exitFullscreen ||
                       document.cancelFullscreen ||  //Why can no one agree on this?
                       document.cancelFullScreen ||
                       document.webkitCancelFullscreen ||
                       document.webkitCancelFullScreen ||
                       document.mozCancelFullscreen ||
                       document.mozCancelFullScreen ||
                       function(){}

Object.defineProperty(proto, "fullscreen", {
  get: function() {
    return this._fullscreenActive
  },
  set: function(state) {
    var ns = !!state
    if(!ns) {
      this._wantFullscreen = false
      cancelFullscreen.call(document)
    } else {
      this._wantFullscreen = true
      tryFullscreen(this)
    }
    return this._fullscreenActive
  }
})

function handleFullscreen(shell) {
  shell._fullscreenActive = document.fullscreen ||
                            document.mozFullScreen ||
                            document.webkitIsFullScreen ||
                            false
  if(!shell.stickyFullscreen && shell._fullscreenActive) {
    shell._wantFullscreen = false
  }
}

//Pointer lock state toggle
var exitPointerLock = document.exitPointerLock ||
                      document.webkitExitPointerLock ||
                      document.mozExitPointerLock ||
                      function() {}

Object.defineProperty(proto, "pointerLock", {
  get: function() {
    return this._pointerLockActive
  },
  set: function(state) {
    var ns = !!state
    if(!ns) {
      this._wantPointerLock = false
      exitPointerLock.call(document)
    } else {
      this._wantPointerLock = true
      tryFullscreen(this)
    }
    return this._pointerLockActive
  }
})

function handlePointerLockChange(shell, event) {
  shell._pointerLockActive = shell.element === (
      document.pointerLockElement ||
      document.mozPointerLockElement ||
      document.webkitPointerLockElement ||
      null)
  if(!shell.stickyPointerLock && shell._pointerLockActive) {
    shell._wantPointerLock = false
  }
}

//Width and height
Object.defineProperty(proto, "width", {
  get: function() {
    return this.element.clientWidth
  }
})
Object.defineProperty(proto, "height", {
  get: function() {
    return this.element.clientHeight
  }
})

//Set key state
function setKeyState(shell, key, state) {
  var ps = shell._curKeyState[key]
  if(ps !== state) {
    if(state) {
      shell._pressCount[key]++
    } else {
      shell._releaseCount[key]++
    }
    shell._curKeyState[key] = state
  }
}

//Ticks the game state one update
function tick(shell) {
  var skip = hrtime() + shell.frameSkip
    , pCount = shell._pressCount
    , rCount = shell._releaseCount
    , i, s, t
    , tr = shell._tickRate
    , n = keyNames.length
  while(!shell._paused &&
        hrtime() >= shell._lastTick + tr) {
    
    //Skip frames if we are over budget
    if(hrtime() > skip) {
      shell._lastTick = hrtime() + tr
      return
    }
    
    //Tick the game
    s = hrtime()
    shell.emit("tick")
    t = hrtime()
    shell.tickTime = t - s
    
    //Update counters and time
    ++shell.tickCount
    shell._lastTick += tr
    
    //Shift input state
    for(i=0; i<n; ++i) {
      pCount[i] = rCount[i] = 0
    }
    if(shell._pointerLockActive) {
      shell.prevMouseX = shell.mouseX = shell.width>>1
      shell.prevMouseY = shell.mouseY = shell.height>>1
    } else {
      shell.prevMouseX = shell.mouseX
      shell.prevMouseY = shell.mouseY
    }
    shell.scroll[0] = shell.scroll[1] = shell.scroll[2] = 0
  }
}

//Render stuff
function render(shell) {

  //Request next frame
  shell._rafHandle = requestAnimationFrame(shell._render)

  //Tick the shell
  tick(shell)
  
  //Compute frame time
  var dt
  if(shell._paused) {
    dt = shell._frameTime
  } else {
    dt = min(1.0, (hrtime() - shell._lastTick) / shell._tickRate)
  }
  
  //Draw a frame
  ++shell.frameCount
  var s = hrtime()
  shell.emit("render", dt)
  var t = hrtime()
  shell.frameTime = t - s
  
}

function isFocused(shell) {
  return (document.activeElement === document.body) ||
         (document.activeElement === shell.element)
}

//Set key up
function handleKeyUp(shell, ev) {
  ev.preventDefault()
  var kc = physicalKeyCode(ev.keyCode || ev.char || ev.which || ev.charCode)
  if(kc >= 0) {
    setKeyState(shell, kc, false)
  }
}

//Set key down
function handleKeyDown(shell, ev) {
  if(!isFocused(shell)) {
    return
  }
  if(ev.metaKey) {
    //Hack: Clear key state when meta gets pressed to prevent keys sticking
    handleBlur(shell, ev)
  } else {
    ev.preventDefault()
    var kc = physicalKeyCode(ev.keyCode || ev.char || ev.which || ev.charCode)
    if(kc >= 0) {
      setKeyState(shell, kc, true)
    }
  }
}

//Mouse events are really annoying
var mouseCodes = iota(32).map(function(n) {
  return virtualKeyCode("mouse-" + (n+1))
})

function setMouseButtons(shell, buttons) {
  for(var i=0; i<32; ++i) {
    setKeyState(shell, mouseCodes[i], !!(buttons & (1<<i)))
  }
}

function handleMouseMove(shell, ev) {
  if(shell._pointerLockActive) {
    var movementX = ev.movementX       ||
                    ev.mozMovementX    ||
                    ev.webkitMovementX ||
                    0,
        movementY = ev.movementY       ||
                    ev.mozMovementY    ||
                    ev.webkitMovementY ||
                    0
    shell.mouseX += movementX
    shell.mouseY += movementY
  } else {
    shell.mouseX = ev.clientX - shell.element.offsetLeft
    shell.mouseY = ev.clientY - shell.element.offsetTop
  }
  return false
}

function handleMouseDown(shell, ev) {
  setKeyState(shell, mouseCodes[ev.button], true)
  return false
}

function handleMouseUp(shell, ev) {
  setKeyState(shell, mouseCodes[ev.button], false)
  return false
}

function handleMouseEnter(shell, ev) {
  if(shell._pointerLockActive) {
    shell.prevMouseX = shell.mouseX = shell.width>>1
    shell.prevMouseY = shell.mouseY = shell.height>>1
  } else {
    shell.prevMouseX = shell.mouseX = ev.clientX - shell.element.offsetLeft
    shell.prevMouseY = shell.mouseY = ev.clientY - shell.element.offsetTop
  }
  return false
}

function handleMouseLeave(shell, ev) {
  setMouseButtons(shell, 0)
  return false
}

//Handle mouse wheel events
function handleMouseWheel(shell, ev) {
  var scale = 1
  switch(ev.deltaMode) {
    case 0: //Pixel
      scale = 1
    break
    case 1: //Line
      scale = 12
    break
    case 2: //Page
       scale = shell.height
    break
  }
  //Add scroll
  shell.scroll[0] +=  ev.deltaX * scale
  shell.scroll[1] +=  ev.deltaY * scale
  shell.scroll[2] += (ev.deltaZ * scale)||0.0
  return false
}

function handleContexMenu(shell, ev) {
  return false
}

function handleBlur(shell, ev) {
  var n = keyNames.length
    , c = shell._curKeyState
    , r = shell._releaseCount
    , i
  for(i=0; i<n; ++i) {
    if(c[i]) {
      ++r[i]
    }
    c[i] = false
  }
  return false
}

function handleResizeElement(shell, ev) {
  var w = shell.element.clientWidth|0
  var h = shell.element.clientHeight|0
  if((w !== shell._width) || (h !== shell._height)) {
    shell._width = w
    shell._height = h
    shell.emit("resize", w, h)
  }
}

function makeDefaultContainer() {
  var container = document.createElement("div")
  container.tabindex = 1
  container.style.position = "absolute"
  container.style.left = "0px"
  container.style.right = "0px"
  container.style.top = "0px"
  container.style.bottom = "0px"
  container.style.height = "100%"
  container.style.overflow = "hidden"
  document.body.appendChild(container)
  document.body.style.overflow = "hidden" //Prevent bounce
  document.body.style.height = "100%"
  return container
}

function createShell(options) {
  options = options || {}
  
  //Check fullscreen and pointer lock flags
  var useFullscreen = !!options.fullscreen
  var usePointerLock = useFullscreen
  if(typeof options.pointerLock !== undefined) {
    usePointerLock = !!options.pointerLock
  }
  
  //Create initial shell
  var shell = new GameShell()
  shell._tickRate = options.tickRate || 30
  shell.frameSkip = options.frameSkip || (shell._tickRate+5) * 5
  shell.stickyFullscreen = !!options.stickyFullscreen || !!options.sticky
  shell.stickyPointerLock = !!options.stickPointerLock || !options.sticky
  
  //Set bindings
  if(options.bindings) {
    shell.bindings = bindings
  }
  
  //Wait for dom to intiailize
  setTimeout(function() { domready(function initGameShell() {
    
    //Retrieve element
    var element = options.element
    if(typeof element === "string") {
      var e = document.querySelector(element)
      if(!e) {
        e = document.getElementById(element)
      }
      if(!e) {
        e = document.getElementByClass(element)[0]
      }
      if(!e) {
        e = makeDefaultContainer()
      }
      shell.element = e
    } else if(typeof element === "object" && !!element) {
      shell.element = element
    } else if(typeof element === "function") {
      shell.element = element()
    } else {
      shell.element = makeDefaultContainer()
    }
    
    //Disable user-select
    if(shell.element.style) {
      shell.element.style["-webkit-touch-callout"] = "none"
      shell.element.style["-webkit-user-select"] = "none"
      shell.element.style["-khtml-user-select"] = "none"
      shell.element.style["-moz-user-select"] = "none"
      shell.element.style["-ms-user-select"] = "none"
      shell.element.style["user-select"] = "none"
    }
    
    //Hook resize handler
    shell._width = shell.element.clientWidth
    shell._height = shell.element.clientHeight
    var handleResize = handleResizeElement.bind(undefined, shell)
    if(typeof MutationObserver !== "undefined") {
      var observer = new MutationObserver(handleResize)
      observer.observe(shell.element, {
        attributes: true,
        subtree: true
      })
    } else {
      shell.element.addEventListener("DOMSubtreeModified", handleResize, false)
    }
    window.addEventListener("resize", handleResize, false)
    
    //Hook keyboard listener
    window.addEventListener("keydown", handleKeyDown.bind(undefined, shell), false)
    window.addEventListener("keyup", handleKeyUp.bind(undefined, shell), false)
    
    //Disable right click
    shell.element.oncontextmenu = handleContexMenu.bind(undefined, shell)
    
    //Hook mouse listeners
    shell.element.addEventListener("mousedown", handleMouseDown.bind(undefined, shell), false)
    shell.element.addEventListener("mouseup", handleMouseUp.bind(undefined, shell), false)
    shell.element.addEventListener("mousemove", handleMouseMove.bind(undefined, shell), false)
    shell.element.addEventListener("mouseenter", handleMouseEnter.bind(undefined, shell), false)
    
    //Mouse leave
    var leave = handleMouseLeave.bind(undefined, shell)
    shell.element.addEventListener("mouseleave", leave, false)
    shell.element.addEventListener("mouseout", leave, false)
    window.addEventListener("mouseleave", leave, false)
    window.addEventListener("mouseout", leave, false)
    
    //Blur event 
    var blur = handleBlur.bind(undefined, shell)
    shell.element.addEventListener("blur", blur, false)
    shell.element.addEventListener("focusout", blur, false)
    shell.element.addEventListener("focus", blur, false)
    window.addEventListener("blur", blur, false)
    window.addEventListener("focusout", blur, false)
    window.addEventListener("focus", blur, false)

    //Mouse wheel handler
    addMouseWheel(shell.element, handleMouseWheel.bind(undefined, shell), false)

    //Fullscreen handler
    var fullscreenChange = handleFullscreen.bind(undefined, shell)
    document.addEventListener("fullscreenchange", fullscreenChange, false)
    document.addEventListener("mozfullscreenchange", fullscreenChange, false)
    document.addEventListener("webkitfullscreenchange", fullscreenChange, false)

    //Stupid fullscreen hack
    shell.element.addEventListener("click", tryFullscreen.bind(undefined, shell), false)

    //Pointer lock change handler
    var pointerLockChange = handlePointerLockChange.bind(undefined, shell)
    document.addEventListener("pointerlockchange", pointerLockChange, false)
    document.addEventListener("mozpointerlockchange", pointerLockChange, false)
    document.addEventListener("webkitpointerlockchange", pointerLockChange, false)
    document.addEventListener("pointerlocklost", pointerLockChange, false)
    document.addEventListener("webkitpointerlocklost", pointerLockChange, false)
    document.addEventListener("mozpointerlocklost", pointerLockChange, false)
    
    //Update flags
    shell.fullscreen = useFullscreen
    shell.pointerLock = usePointerLock
  
    //Default mouse button aliases
    shell.bind("mouse-left",   "mouse-1")
    shell.bind("mouse-right",  "mouse-3")
    shell.bind("mouse-middle", "mouse-2")
    
    //Initialize tick counter
    shell._lastTick = hrtime()
    shell.startTime = hrtime()

    //Unpause shell
    shell.paused = false
    
    //Emit initialize event
    shell.emit("init")
  })}, 0)
  
  return shell
}

module.exports = createShell

},{"./lib/hrtime-polyfill.js":"/Users/vallette/Desktop/colony/node_modules/game-shell/lib/hrtime-polyfill.js","./lib/mousewheel-polyfill.js":"/Users/vallette/Desktop/colony/node_modules/game-shell/lib/mousewheel-polyfill.js","./lib/raf-polyfill.js":"/Users/vallette/Desktop/colony/node_modules/game-shell/lib/raf-polyfill.js","binary-search-bounds":"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/binary-search-bounds/search-bounds.js","domready":"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/domready/ready.js","events":"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/events/events.js","invert-hash":"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/invert-hash/invert.js","iota-array":"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/iota-array/iota.js","uniq":"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/uniq/uniq.js","util":"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/util/util.js","vkey":"/Users/vallette/Desktop/colony/node_modules/game-shell/node_modules/vkey/index.js"}],"/Users/vallette/Desktop/colony/node_modules/harmony-collections/harmony-collections.js":[function(require,module,exports){
/* (The MIT License)
 *
 * Copyright (c) 2012 Brandon Benvie <http://bbenvie.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the 'Software'), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included with all copies or
 * substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 * BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

// Original WeakMap implementation by Gozala @ https://gist.github.com/1269991
// Updated and bugfixed by Raynos @ https://gist.github.com/1638059
// Expanded by Benvie @ https://github.com/Benvie/harmony-collections

void function(string_, object_, function_, prototype_, toString_,
              Array, Object, Function, FP, global, exports, undefined_, undefined){

  var getProperties = Object.getOwnPropertyNames,
      es5 = typeof getProperties === function_ && !(prototype_ in getProperties);

  var callbind = FP.bind
    ? FP.bind.bind(FP.call)
    : (function(call){
        return function(func){
          return function(){
            return call.apply(func, arguments);
          };
        };
      }(FP.call));

  var functionToString = callbind(FP[toString_]),
      objectToString = callbind({}[toString_]),
      numberToString = callbind(.0.toString),
      call = callbind(FP.call),
      apply = callbind(FP.apply),
      hasOwn = callbind({}.hasOwnProperty),
      push = callbind([].push),
      splice = callbind([].splice);

  var name = function(func){
    if (typeof func !== function_)
      return '';
    else if ('name' in func)
      return func.name;

    return functionToString(func).match(/^\n?function\s?(\w*)?_?\(/)[1];
  };

  var create = es5
    ? Object.create
    : function(proto, descs){
        var Ctor = function(){};
        Ctor[prototype_] = Object(proto);
        var object = new Ctor;

        if (descs)
          for (var key in descs)
            defineProperty(object, key, descs[k]);

        return object;
      };


  function Hash(){}

  if (es5) {
    void function(ObjectCreate){
      Hash.prototype = ObjectCreate(null);
      function inherit(obj){
        return ObjectCreate(obj);
      }
      Hash.inherit = inherit;
    }(Object.create);
  } else {
    void function(F){
      var iframe = document.createElement('iframe');
      iframe.style.display = 'none';
      document.body.appendChild(iframe);
      iframe.src = 'javascript:'
      Hash.prototype = iframe.contentWindow.Object.prototype;
      document.body.removeChild(iframe);
      iframe = null;

      var props = ['constructor', 'hasOwnProperty', 'propertyIsEnumerable',
                   'isProtoypeOf', 'toLocaleString', 'toString', 'valueOf'];

      for (var i=0; i < props.length; i++)
        delete Hash.prototype[props[i]];

      function inherit(obj){
        F.prototype = obj;
        obj = new F;
        F.prototype = null;
        return obj;
      }

      Hash.inherit = inherit;
    }(function(){});
  }

  var defineProperty = es5
    ? Object.defineProperty
    : function(object, key, desc) {
        object[key] = desc.value;
        return object;
      };

  var define = function(object, key, value){
    if (typeof key === function_) {
      value = key;
      key = name(value).replace(/_$/, '');
    }

    return defineProperty(object, key, { configurable: true, writable: true, value: value });
  };

  var isArray = es5
    ? (function(isArray){
        return function(o){
          return isArray(o) || o instanceof Array;
        };
      })(Array.isArray)
    : function(o){
        return o instanceof Array || objectToString(o) === '[object Array]';
      };

  // ############
  // ### Data ###
  // ############

  var builtinWeakMap = 'WeakMap' in global;

  var MapData = builtinWeakMap
    ? (function(){
      var BuiltinWeakMap = global.WeakMap,
          wmget = callbind(BuiltinWeakMap[prototype_].get),
          wmset = callbind(BuiltinWeakMap[prototype_].set),
          wmhas = callbind(BuiltinWeakMap[prototype_].has);

      function MapData(name){
        var map = new BuiltinWeakMap;

        this.get = function(o){
          return wmget(map, o);
        };
        this.set = function(o, v){
          wmset(map, o, v);
        };

        if (name) {
          this.wrap = function(o, v){
            if (wmhas(map, o))
              throw new TypeError("Object is already a " + name);
            wmset(map, o, v);
          };
          this.unwrap = function(o){
            var storage = wmget(map, o);
            if (!storage)
              throw new TypeError(name + " is not generic");
            return storage;
          };
        }
      }

      return MapData;
    })()
    : (function(){
      var locker = 'return function(k){if(k===s)return l}',
          random = Math.random,
          uids = new Hash,
          slice = callbind(''.slice),
          indexOf = callbind([].indexOf);

      var createUID = function(){
        var key = slice(numberToString(random(), 36), 2);
        return key in uids ? createUID() : uids[key] = key;
      };

      var globalID = createUID();

      // common per-object storage area made visible by patching getOwnPropertyNames'
      function getOwnPropertyNames(obj){
        var props = getProperties(obj);
        if (hasOwn(obj, globalID))
          splice(props, indexOf(props, globalID), 1);
        return props;
      }

      if (es5) {
        // check for the random key on an object, create new storage if missing, return it
        var storage = function(obj){
          if (!hasOwn(obj, globalID))
            defineProperty(obj, globalID, { value: new Hash });
          return obj[globalID];
        };

        define(Object, getOwnPropertyNames);
      } else {

        var toStringToString = function(s){
          function toString(){ return s }
          return toString[toString_] = toString;
        }(Object[prototype_][toString_]+'');

        // store the values on a custom valueOf in order to hide them but store them locally
        var storage = function(obj){
          if (hasOwn(obj, toString_) && globalID in obj[toString_])
            return obj[toString_][globalID];

          if (!(toString_ in obj))
            throw new Error("Can't store values for "+obj);

          var oldToString = obj[toString_];
          function toString(){ return oldToString.call(this) }
          obj[toString_] = toString;
          toString[toString_] = toStringToString;
          return toString[globalID] = {};
        };
      }



      // shim for [[MapData]] from es6 spec, and pulls double duty as WeakMap storage
      function MapData(name){
        var puid = createUID(),
            iuid = createUID(),
            secret = { value: undefined };

        var attach = function(obj){
          var store = storage(obj);
          if (hasOwn(store, puid))
            return store[puid](secret);

          var lockbox = new Hash;
          defineProperty(lockbox, iuid, secret);
          defineProperty(store, puid, {
            value: new Function('s', 'l', locker)(secret, lockbox)
          });
          return lockbox;
        };

        this.get = function(o){
          return attach(o)[iuid];
        };
        this.set = function(o, v){
          attach(o)[iuid] = v;
        };

        if (name) {
          this.wrap = function(o, v){
            var lockbox = attach(o);
            if (lockbox[iuid])
              throw new TypeError("Object is already a " + name);
            lockbox[iuid] = v;
          };
          this.unwrap = function(o){
            var storage = attach(o)[iuid];
            if (!storage)
              throw new TypeError(name + " is not generic");
            return storage;
          };
        }
      }

      return MapData;
    }());

  var exporter = (function(){
    // [native code] looks slightly different in each engine
    var src = (''+Object).split('Object');

    // fake [native code]
    function toString(){
      return src[0] + name(this) + src[1];
    }

    define(toString, toString);

    // attempt to use __proto__ so the methods don't all have an own toString
    var prepFunction = { __proto__: [] } instanceof Array
      ? function(func){ func.__proto__ = toString }
      : function(func){ define(func, toString) };

    // assemble an array of functions into a fully formed class
    var prepare = function(methods){
      var Ctor = methods.shift(),
          brand = '[object ' + name(Ctor) + ']';

      function toString(){ return brand }
      methods.push(toString);
      prepFunction(Ctor);

      for (var i=0; i < methods.length; i++) {
        prepFunction(methods[i]);
        define(Ctor[prototype_], methods[i]);
      }

      return Ctor;
    };

    return function(name, init){
      if (name in exports)
        return exports[name];

      var data = new MapData(name);

      return exports[name] = prepare(init(
        function(collection, value){
          data.wrap(collection, value);
        },
        function(collection){
          return data.unwrap(collection);
        }
      ));
    };
  }());


  // initialize collection with an iterable, currently only supports forEach function
  var initialize = function(iterable, callback){
    if (iterable !== null && typeof iterable === object_ && typeof iterable.forEach === function_) {
      iterable.forEach(function(item, i){
        if (isArray(item) && item.length === 2)
          callback(iterable[i][0], iterable[i][1]);
        else
          callback(iterable[i], i);
      });
    }
  }

  // attempt to fix the name of "delete_" methods, should work in V8 and spidermonkey
  var fixDelete = function(func, scopeNames, scopeValues){
    try {
      scopeNames[scopeNames.length] = ('return '+func).replace('e_', '\\u0065');
      return Function.apply(0, scopeNames).apply(0, scopeValues);
    } catch (e) {
      return func;
    }
  }

  var WM, HM, M;

  // ###############
  // ### WeakMap ###
  // ###############

  WM = builtinWeakMap ? (exports.WeakMap = global.WeakMap) : exporter('WeakMap', function(wrap, unwrap){
    var prototype = WeakMap[prototype_];
    var validate = function(key){
      if (key == null || typeof key !== object_ && typeof key !== function_)
        throw new TypeError("Invalid WeakMap key");
    };

    /**
     * @class        WeakMap
     * @description  Collection using objects with unique identities as keys that disallows enumeration
     *               and allows for better garbage collection.
     * @param        {Iterable} [iterable]  An item to populate the collection with.
     */
    function WeakMap(iterable){
      if (this === global || this == null || this === prototype)
        return new WeakMap(iterable);

      wrap(this, new MapData);

      var self = this;
      iterable && initialize(iterable, function(value, key){
        call(set, self, value, key);
      });
    }
    /**
     * @method       <get>
     * @description  Retrieve the value in the collection that matches key
     * @param        {Any} key
     * @return       {Any}
     */
    function get(key){
      validate(key);
      var value = unwrap(this).get(key);
      return value === undefined_ ? undefined : value;
    }
    /**
     * @method       <set>
     * @description  Add or update a pair in the collection. Enforces uniqueness by overwriting.
     * @param        {Any} key
     * @param        {Any} val
     **/
    function set(key, value){
      validate(key);
      // store a token for explicit undefined so that "has" works correctly
      unwrap(this).set(key, value === undefined ? undefined_ : value);
    }
    /*
     * @method       <has>
     * @description  Check if key is in the collection
     * @param        {Any} key
     * @return       {Boolean}
     **/
    function has(key){
      validate(key);
      return unwrap(this).get(key) !== undefined;
    }
    /**
     * @method       <delete>
     * @description  Remove key and matching value if found
     * @param        {Any} key
     * @return       {Boolean} true if item was in collection
     */
    function delete_(key){
      validate(key);
      var data = unwrap(this);

      if (data.get(key) === undefined)
        return false;

      data.set(key, undefined);
      return true;
    }

    delete_ = fixDelete(delete_, ['validate', 'unwrap'], [validate, unwrap]);
    return [WeakMap, get, set, has, delete_];
  });


  // ###############
  // ### HashMap ###
  // ###############

  HM = exporter('HashMap', function(wrap, unwrap){
    // separate numbers, strings, and atoms to compensate for key coercion to string

    var prototype = HashMap[prototype_],
        STRING = 0, NUMBER = 1, OTHER = 2,
        others = { 'true': true, 'false': false, 'null': null, 0: -0 };

    var proto = Math.random().toString(36).slice(2);

    var coerce = function(key){
      return key === '__proto__' ? proto : key;
    };

    var uncoerce = function(type, key){
      switch (type) {
        case STRING: return key === proto ? '__proto__' : key;
        case NUMBER: return +key;
        case OTHER: return others[key];
      }
    }


    var validate = function(key){
      if (key == null) return OTHER;
      switch (typeof key) {
        case 'boolean': return OTHER;
        case string_: return STRING;
        // negative zero has to be explicitly accounted for
        case 'number': return key === 0 && Infinity / key === -Infinity ? OTHER : NUMBER;
        default: throw new TypeError("Invalid HashMap key");
      }
    }

    /**
     * @class          HashMap
     * @description    Collection that only allows primitives to be keys.
     * @param          {Iterable} [iterable]  An item to populate the collection with.
     */
    function HashMap(iterable){
      if (this === global || this == null || this === prototype)
        return new HashMap(iterable);

      wrap(this, {
        size: 0,
        0: new Hash,
        1: new Hash,
        2: new Hash
      });

      var self = this;
      iterable && initialize(iterable, function(value, key){
        call(set, self, value, key);
      });
    }
    /**
     * @method       <get>
     * @description  Retrieve the value in the collection that matches key
     * @param        {Any} key
     * @return       {Any}
     */
    function get(key){
      return unwrap(this)[validate(key)][coerce(key)];
    }
    /**
     * @method       <set>
     * @description  Add or update a pair in the collection. Enforces uniqueness by overwriting.
     * @param        {Any} key
     * @param        {Any} val
     **/
    function set(key, value){
      var items = unwrap(this),
          data = items[validate(key)];

      key = coerce(key);
      key in data || items.size++;
      data[key] = value;
    }
    /**
     * @method       <has>
     * @description  Check if key exists in the collection.
     * @param        {Any} key
     * @return       {Boolean} is in collection
     **/
    function has(key){
      return coerce(key) in unwrap(this)[validate(key)];
    }
    /**
     * @method       <delete>
     * @description  Remove key and matching value if found
     * @param        {Any} key
     * @return       {Boolean} true if item was in collection
     */
    function delete_(key){
      var items = unwrap(this),
          data = items[validate(key)];

      key = coerce(key);
      if (key in data) {
        delete data[key];
        items.size--;
        return true;
      }

      return false;
    }
    /**
     * @method       <size>
     * @description  Retrieve the amount of items in the collection
     * @return       {Number}
     */
    function size(){
      return unwrap(this).size;
    }
    /**
     * @method       <forEach>
     * @description  Loop through the collection raising callback for each
     * @param        {Function} callback  `callback(value, key)`
     * @param        {Object}   context    The `this` binding for callbacks, default null
     */
    function forEach(callback, context){
      var data = unwrap(this);
      context = context == null ? global : context;
      for (var i=0; i < 3; i++)
        for (var key in data[i])
          call(callback, context, data[i][key], uncoerce(i, key), this);
    }

    delete_ = fixDelete(delete_, ['validate', 'unwrap', 'coerce'], [validate, unwrap, coerce]);
    return [HashMap, get, set, has, delete_, size, forEach];
  });


  // ###########
  // ### Map ###
  // ###########

  // if a fully implemented Map exists then use it
  if ('Map' in global && 'forEach' in global.Map.prototype) {
    M = exports.Map = global.Map;
  } else {
    M = exporter('Map', function(wrap, unwrap){
      // attempt to use an existing partially implemented Map
      var BuiltinMap = global.Map,
          prototype = Map[prototype_],
          wm = WM[prototype_],
          hm = (BuiltinMap || HM)[prototype_],
          mget    = [callbind(hm.get), callbind(wm.get)],
          mset    = [callbind(hm.set), callbind(wm.set)],
          mhas    = [callbind(hm.has), callbind(wm.has)],
          mdelete = [callbind(hm['delete']), callbind(wm['delete'])];

      var type = BuiltinMap
        ? function(){ return 0 }
        : function(o){ return +(typeof o === object_ ? o !== null : typeof o === function_) }

      // if we have a builtin Map we can let it do most of the heavy lifting
      var init = BuiltinMap
        ? function(){ return { 0: new BuiltinMap } }
        : function(){ return { 0: new HM, 1: new WM } };

      /**
       * @class         Map
       * @description   Collection that allows any kind of value to be a key.
       * @param         {Iterable} [iterable]  An item to populate the collection with.
       */
      function Map(iterable){
        if (this === global || this == null || this === prototype)
          return new Map(iterable);

        var data = init();
        data.keys = [];
        data.values = [];
        wrap(this, data);

        var self = this;
        iterable && initialize(iterable, function(value, key){
          call(set, self, value, key);
        });
      }
      /**
       * @method       <get>
       * @description  Retrieve the value in the collection that matches key
       * @param        {Any} key
       * @return       {Any}
       */
      function get(key){
        var data = unwrap(this),
            t = type(key);
        return data.values[mget[t](data[t], key)];
      }
      /**
       * @method       <set>
       * @description  Add or update a pair in the collection. Enforces uniqueness by overwriting.
       * @param        {Any} key
       * @param        {Any} val
       **/
      function set(key, value){
        var data = unwrap(this),
            t = type(key),
            index = mget[t](data[t], key);

        if (index === undefined) {
          mset[t](data[t], key, data.keys.length);
          push(data.keys, key);
          push(data.values, value);
        } else {
          data.keys[index] = key;
          data.values[index] = value;
        }
      }
      /**
       * @method       <has>
       * @description  Check if key exists in the collection.
       * @param        {Any} key
       * @return       {Boolean} is in collection
       **/
      function has(key){
        var t = type(key);
        return mhas[t](unwrap(this)[t], key);
      }
      /**
       * @method       <delete>
       * @description  Remove key and matching value if found
       * @param        {Any} key
       * @return       {Boolean} true if item was in collection
       */
      function delete_(key){
        var data = unwrap(this),
            t = type(key),
            index = mget[t](data[t], key);

        if (index === undefined)
          return false;

        mdelete[t](data[t], key);
        splice(data.keys, index, 1);
        splice(data.values, index, 1);
        return true;
      }
      /**
       * @method       <size>
       * @description  Retrieve the amount of items in the collection
       * @return       {Number}
       */
      function size(){
        return unwrap(this).keys.length;
      }
      /**
       * @method       <forEach>
       * @description  Loop through the collection raising callback for each
       * @param        {Function} callback  `callback(value, key)`
       * @param        {Object}   context    The `this` binding for callbacks, default null
       */
      function forEach(callback, context){
        var data = unwrap(this),
            keys = data.keys,
            values = data.values;

        context = context == null ? global : context;

        for (var i=0, len=keys.length; i < len; i++)
          call(callback, context, values[i], keys[i], this);
      }

      delete_ = fixDelete(delete_,
        ['type', 'unwrap', 'call', 'splice'],
        [type, unwrap, call, splice]
      );
      return [Map, get, set, has, delete_, size, forEach];
    });
  }


  // ###########
  // ### Set ###
  // ###########

  exporter('Set', function(wrap, unwrap){
    var prototype = Set[prototype_],
        m = M[prototype_],
        msize = callbind(m.size),
        mforEach = callbind(m.forEach),
        mget = callbind(m.get),
        mset = callbind(m.set),
        mhas = callbind(m.has),
        mdelete = callbind(m['delete']);

    /**
     * @class        Set
     * @description  Collection of values that enforces uniqueness.
     * @param        {Iterable} [iterable]  An item to populate the collection with.
     **/
    function Set(iterable){
      if (this === global || this == null || this === prototype)
        return new Set(iterable);

      wrap(this, new M);

      var self = this;
      iterable && initialize(iterable, function(value, key){
        call(add, self, key);
      });
    }
    /**
     * @method       <add>
     * @description  Insert value if not found, enforcing uniqueness.
     * @param        {Any} val
     */
    function add(key){
      mset(unwrap(this), key, key);
    }
    /**
     * @method       <has>
     * @description  Check if key exists in the collection.
     * @param        {Any} key
     * @return       {Boolean} is in collection
     **/
    function has(key){
      return mhas(unwrap(this), key);
    }
    /**
     * @method       <delete>
     * @description  Remove key and matching value if found
     * @param        {Any} key
     * @return       {Boolean} true if item was in collection
     */
    function delete_(key){
      return mdelete(unwrap(this), key);
    }
    /**
     * @method       <size>
     * @description  Retrieve the amount of items in the collection
     * @return       {Number}
     */
    function size(){
      return msize(unwrap(this));
    }
    /**
     * @method       <forEach>
     * @description  Loop through the collection raising callback for each. Index is simply the counter for the current iteration.
     * @param        {Function} callback  `callback(value, index)`
     * @param        {Object}   context    The `this` binding for callbacks, default null
     */
    function forEach(callback, context){
      var index = 0,
          self = this;
      mforEach(unwrap(this, function(key){
        call(callback, this, key, index++, self);
      }, context));
    }

    delete_ = fixDelete(delete_, ['mdelete', 'unwrap'], [mdelete, unwrap]);
    return [Set, add, has, delete_, size, forEach];
  });
}('string', 'object', 'function', 'prototype', 'toString',
  Array, Object, Function, Function.prototype, (0, eval)('this'),
  typeof exports === 'undefined' ? this : exports, {});

},{}],"/Users/vallette/Desktop/colony/src/index.js":[function(require,module,exports){
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
},{"delaunay-triangulate":"/Users/vallette/Desktop/colony/node_modules/delaunay-triangulate/triangulate.js","game-shell":"/Users/vallette/Desktop/colony/node_modules/game-shell/shell.js","harmony-collections":"/Users/vallette/Desktop/colony/node_modules/harmony-collections/harmony-collections.js"}],"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/events/events.js":[function(require,module,exports){
// Copyright Joyent, Inc. and other Node contributors.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit
// persons to whom the Software is furnished to do so, subject to the
// following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
// NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
// USE OR OTHER DEALINGS IN THE SOFTWARE.

function EventEmitter() {
  this._events = this._events || {};
  this._maxListeners = this._maxListeners || undefined;
}
module.exports = EventEmitter;

// Backwards-compat with node 0.10.x
EventEmitter.EventEmitter = EventEmitter;

EventEmitter.prototype._events = undefined;
EventEmitter.prototype._maxListeners = undefined;

// By default EventEmitters will print a warning if more than 10 listeners are
// added to it. This is a useful default which helps finding memory leaks.
EventEmitter.defaultMaxListeners = 10;

// Obviously not all Emitters should be limited to 10. This function allows
// that to be increased. Set to zero for unlimited.
EventEmitter.prototype.setMaxListeners = function(n) {
  if (!isNumber(n) || n < 0 || isNaN(n))
    throw TypeError('n must be a positive number');
  this._maxListeners = n;
  return this;
};

EventEmitter.prototype.emit = function(type) {
  var er, handler, len, args, i, listeners;

  if (!this._events)
    this._events = {};

  // If there is no 'error' event listener then throw.
  if (type === 'error') {
    if (!this._events.error ||
        (isObject(this._events.error) && !this._events.error.length)) {
      er = arguments[1];
      if (er instanceof Error) {
        throw er; // Unhandled 'error' event
      } else {
        throw TypeError('Uncaught, unspecified "error" event.');
      }
      return false;
    }
  }

  handler = this._events[type];

  if (isUndefined(handler))
    return false;

  if (isFunction(handler)) {
    switch (arguments.length) {
      // fast cases
      case 1:
        handler.call(this);
        break;
      case 2:
        handler.call(this, arguments[1]);
        break;
      case 3:
        handler.call(this, arguments[1], arguments[2]);
        break;
      // slower
      default:
        len = arguments.length;
        args = new Array(len - 1);
        for (i = 1; i < len; i++)
          args[i - 1] = arguments[i];
        handler.apply(this, args);
    }
  } else if (isObject(handler)) {
    len = arguments.length;
    args = new Array(len - 1);
    for (i = 1; i < len; i++)
      args[i - 1] = arguments[i];

    listeners = handler.slice();
    len = listeners.length;
    for (i = 0; i < len; i++)
      listeners[i].apply(this, args);
  }

  return true;
};

EventEmitter.prototype.addListener = function(type, listener) {
  var m;

  if (!isFunction(listener))
    throw TypeError('listener must be a function');

  if (!this._events)
    this._events = {};

  // To avoid recursion in the case that type === "newListener"! Before
  // adding it to the listeners, first emit "newListener".
  if (this._events.newListener)
    this.emit('newListener', type,
              isFunction(listener.listener) ?
              listener.listener : listener);

  if (!this._events[type])
    // Optimize the case of one listener. Don't need the extra array object.
    this._events[type] = listener;
  else if (isObject(this._events[type]))
    // If we've already got an array, just append.
    this._events[type].push(listener);
  else
    // Adding the second element, need to change to array.
    this._events[type] = [this._events[type], listener];

  // Check for listener leak
  if (isObject(this._events[type]) && !this._events[type].warned) {
    var m;
    if (!isUndefined(this._maxListeners)) {
      m = this._maxListeners;
    } else {
      m = EventEmitter.defaultMaxListeners;
    }

    if (m && m > 0 && this._events[type].length > m) {
      this._events[type].warned = true;
      console.error('(node) warning: possible EventEmitter memory ' +
                    'leak detected. %d listeners added. ' +
                    'Use emitter.setMaxListeners() to increase limit.',
                    this._events[type].length);
      if (typeof console.trace === 'function') {
        // not supported in IE 10
        console.trace();
      }
    }
  }

  return this;
};

EventEmitter.prototype.on = EventEmitter.prototype.addListener;

EventEmitter.prototype.once = function(type, listener) {
  if (!isFunction(listener))
    throw TypeError('listener must be a function');

  var fired = false;

  function g() {
    this.removeListener(type, g);

    if (!fired) {
      fired = true;
      listener.apply(this, arguments);
    }
  }

  g.listener = listener;
  this.on(type, g);

  return this;
};

// emits a 'removeListener' event iff the listener was removed
EventEmitter.prototype.removeListener = function(type, listener) {
  var list, position, length, i;

  if (!isFunction(listener))
    throw TypeError('listener must be a function');

  if (!this._events || !this._events[type])
    return this;

  list = this._events[type];
  length = list.length;
  position = -1;

  if (list === listener ||
      (isFunction(list.listener) && list.listener === listener)) {
    delete this._events[type];
    if (this._events.removeListener)
      this.emit('removeListener', type, listener);

  } else if (isObject(list)) {
    for (i = length; i-- > 0;) {
      if (list[i] === listener ||
          (list[i].listener && list[i].listener === listener)) {
        position = i;
        break;
      }
    }

    if (position < 0)
      return this;

    if (list.length === 1) {
      list.length = 0;
      delete this._events[type];
    } else {
      list.splice(position, 1);
    }

    if (this._events.removeListener)
      this.emit('removeListener', type, listener);
  }

  return this;
};

EventEmitter.prototype.removeAllListeners = function(type) {
  var key, listeners;

  if (!this._events)
    return this;

  // not listening for removeListener, no need to emit
  if (!this._events.removeListener) {
    if (arguments.length === 0)
      this._events = {};
    else if (this._events[type])
      delete this._events[type];
    return this;
  }

  // emit removeListener for all listeners on all events
  if (arguments.length === 0) {
    for (key in this._events) {
      if (key === 'removeListener') continue;
      this.removeAllListeners(key);
    }
    this.removeAllListeners('removeListener');
    this._events = {};
    return this;
  }

  listeners = this._events[type];

  if (isFunction(listeners)) {
    this.removeListener(type, listeners);
  } else {
    // LIFO order
    while (listeners.length)
      this.removeListener(type, listeners[listeners.length - 1]);
  }
  delete this._events[type];

  return this;
};

EventEmitter.prototype.listeners = function(type) {
  var ret;
  if (!this._events || !this._events[type])
    ret = [];
  else if (isFunction(this._events[type]))
    ret = [this._events[type]];
  else
    ret = this._events[type].slice();
  return ret;
};

EventEmitter.listenerCount = function(emitter, type) {
  var ret;
  if (!emitter._events || !emitter._events[type])
    ret = 0;
  else if (isFunction(emitter._events[type]))
    ret = 1;
  else
    ret = emitter._events[type].length;
  return ret;
};

function isFunction(arg) {
  return typeof arg === 'function';
}

function isNumber(arg) {
  return typeof arg === 'number';
}

function isObject(arg) {
  return typeof arg === 'object' && arg !== null;
}

function isUndefined(arg) {
  return arg === void 0;
}

},{}],"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/inherits/inherits_browser.js":[function(require,module,exports){
if (typeof Object.create === 'function') {
  // implementation from standard node.js 'util' module
  module.exports = function inherits(ctor, superCtor) {
    ctor.super_ = superCtor
    ctor.prototype = Object.create(superCtor.prototype, {
      constructor: {
        value: ctor,
        enumerable: false,
        writable: true,
        configurable: true
      }
    });
  };
} else {
  // old school shim for old browsers
  module.exports = function inherits(ctor, superCtor) {
    ctor.super_ = superCtor
    var TempCtor = function () {}
    TempCtor.prototype = superCtor.prototype
    ctor.prototype = new TempCtor()
    ctor.prototype.constructor = ctor
  }
}

},{}],"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/process/browser.js":[function(require,module,exports){
// shim for using process in browser

var process = module.exports = {};

process.nextTick = (function () {
    var canSetImmediate = typeof window !== 'undefined'
    && window.setImmediate;
    var canPost = typeof window !== 'undefined'
    && window.postMessage && window.addEventListener
    ;

    if (canSetImmediate) {
        return function (f) { return window.setImmediate(f) };
    }

    if (canPost) {
        var queue = [];
        window.addEventListener('message', function (ev) {
            var source = ev.source;
            if ((source === window || source === null) && ev.data === 'process-tick') {
                ev.stopPropagation();
                if (queue.length > 0) {
                    var fn = queue.shift();
                    fn();
                }
            }
        }, true);

        return function nextTick(fn) {
            queue.push(fn);
            window.postMessage('process-tick', '*');
        };
    }

    return function nextTick(fn) {
        setTimeout(fn, 0);
    };
})();

process.title = 'browser';
process.browser = true;
process.env = {};
process.argv = [];

function noop() {}

process.on = noop;
process.addListener = noop;
process.once = noop;
process.off = noop;
process.removeListener = noop;
process.removeAllListeners = noop;
process.emit = noop;

process.binding = function (name) {
    throw new Error('process.binding is not supported');
}

// TODO(shtylman)
process.cwd = function () { return '/' };
process.chdir = function (dir) {
    throw new Error('process.chdir is not supported');
};

},{}],"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/util/support/isBufferBrowser.js":[function(require,module,exports){
module.exports = function isBuffer(arg) {
  return arg && typeof arg === 'object'
    && typeof arg.copy === 'function'
    && typeof arg.fill === 'function'
    && typeof arg.readUInt8 === 'function';
}
},{}],"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/util/util.js":[function(require,module,exports){
(function (process,global){
// Copyright Joyent, Inc. and other Node contributors.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit
// persons to whom the Software is furnished to do so, subject to the
// following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
// NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
// USE OR OTHER DEALINGS IN THE SOFTWARE.

var formatRegExp = /%[sdj%]/g;
exports.format = function(f) {
  if (!isString(f)) {
    var objects = [];
    for (var i = 0; i < arguments.length; i++) {
      objects.push(inspect(arguments[i]));
    }
    return objects.join(' ');
  }

  var i = 1;
  var args = arguments;
  var len = args.length;
  var str = String(f).replace(formatRegExp, function(x) {
    if (x === '%%') return '%';
    if (i >= len) return x;
    switch (x) {
      case '%s': return String(args[i++]);
      case '%d': return Number(args[i++]);
      case '%j':
        try {
          return JSON.stringify(args[i++]);
        } catch (_) {
          return '[Circular]';
        }
      default:
        return x;
    }
  });
  for (var x = args[i]; i < len; x = args[++i]) {
    if (isNull(x) || !isObject(x)) {
      str += ' ' + x;
    } else {
      str += ' ' + inspect(x);
    }
  }
  return str;
};


// Mark that a method should not be used.
// Returns a modified function which warns once by default.
// If --no-deprecation is set, then it is a no-op.
exports.deprecate = function(fn, msg) {
  // Allow for deprecating things in the process of starting up.
  if (isUndefined(global.process)) {
    return function() {
      return exports.deprecate(fn, msg).apply(this, arguments);
    };
  }

  if (process.noDeprecation === true) {
    return fn;
  }

  var warned = false;
  function deprecated() {
    if (!warned) {
      if (process.throwDeprecation) {
        throw new Error(msg);
      } else if (process.traceDeprecation) {
        console.trace(msg);
      } else {
        console.error(msg);
      }
      warned = true;
    }
    return fn.apply(this, arguments);
  }

  return deprecated;
};


var debugs = {};
var debugEnviron;
exports.debuglog = function(set) {
  if (isUndefined(debugEnviron))
    debugEnviron = process.env.NODE_DEBUG || '';
  set = set.toUpperCase();
  if (!debugs[set]) {
    if (new RegExp('\\b' + set + '\\b', 'i').test(debugEnviron)) {
      var pid = process.pid;
      debugs[set] = function() {
        var msg = exports.format.apply(exports, arguments);
        console.error('%s %d: %s', set, pid, msg);
      };
    } else {
      debugs[set] = function() {};
    }
  }
  return debugs[set];
};


/**
 * Echos the value of a value. Trys to print the value out
 * in the best way possible given the different types.
 *
 * @param {Object} obj The object to print out.
 * @param {Object} opts Optional options object that alters the output.
 */
/* legacy: obj, showHidden, depth, colors*/
function inspect(obj, opts) {
  // default options
  var ctx = {
    seen: [],
    stylize: stylizeNoColor
  };
  // legacy...
  if (arguments.length >= 3) ctx.depth = arguments[2];
  if (arguments.length >= 4) ctx.colors = arguments[3];
  if (isBoolean(opts)) {
    // legacy...
    ctx.showHidden = opts;
  } else if (opts) {
    // got an "options" object
    exports._extend(ctx, opts);
  }
  // set default options
  if (isUndefined(ctx.showHidden)) ctx.showHidden = false;
  if (isUndefined(ctx.depth)) ctx.depth = 2;
  if (isUndefined(ctx.colors)) ctx.colors = false;
  if (isUndefined(ctx.customInspect)) ctx.customInspect = true;
  if (ctx.colors) ctx.stylize = stylizeWithColor;
  return formatValue(ctx, obj, ctx.depth);
}
exports.inspect = inspect;


// http://en.wikipedia.org/wiki/ANSI_escape_code#graphics
inspect.colors = {
  'bold' : [1, 22],
  'italic' : [3, 23],
  'underline' : [4, 24],
  'inverse' : [7, 27],
  'white' : [37, 39],
  'grey' : [90, 39],
  'black' : [30, 39],
  'blue' : [34, 39],
  'cyan' : [36, 39],
  'green' : [32, 39],
  'magenta' : [35, 39],
  'red' : [31, 39],
  'yellow' : [33, 39]
};

// Don't use 'blue' not visible on cmd.exe
inspect.styles = {
  'special': 'cyan',
  'number': 'yellow',
  'boolean': 'yellow',
  'undefined': 'grey',
  'null': 'bold',
  'string': 'green',
  'date': 'magenta',
  // "name": intentionally not styling
  'regexp': 'red'
};


function stylizeWithColor(str, styleType) {
  var style = inspect.styles[styleType];

  if (style) {
    return '\u001b[' + inspect.colors[style][0] + 'm' + str +
           '\u001b[' + inspect.colors[style][1] + 'm';
  } else {
    return str;
  }
}


function stylizeNoColor(str, styleType) {
  return str;
}


function arrayToHash(array) {
  var hash = {};

  array.forEach(function(val, idx) {
    hash[val] = true;
  });

  return hash;
}


function formatValue(ctx, value, recurseTimes) {
  // Provide a hook for user-specified inspect functions.
  // Check that value is an object with an inspect function on it
  if (ctx.customInspect &&
      value &&
      isFunction(value.inspect) &&
      // Filter out the util module, it's inspect function is special
      value.inspect !== exports.inspect &&
      // Also filter out any prototype objects using the circular check.
      !(value.constructor && value.constructor.prototype === value)) {
    var ret = value.inspect(recurseTimes, ctx);
    if (!isString(ret)) {
      ret = formatValue(ctx, ret, recurseTimes);
    }
    return ret;
  }

  // Primitive types cannot have properties
  var primitive = formatPrimitive(ctx, value);
  if (primitive) {
    return primitive;
  }

  // Look up the keys of the object.
  var keys = Object.keys(value);
  var visibleKeys = arrayToHash(keys);

  if (ctx.showHidden) {
    keys = Object.getOwnPropertyNames(value);
  }

  // IE doesn't make error fields non-enumerable
  // http://msdn.microsoft.com/en-us/library/ie/dww52sbt(v=vs.94).aspx
  if (isError(value)
      && (keys.indexOf('message') >= 0 || keys.indexOf('description') >= 0)) {
    return formatError(value);
  }

  // Some type of object without properties can be shortcutted.
  if (keys.length === 0) {
    if (isFunction(value)) {
      var name = value.name ? ': ' + value.name : '';
      return ctx.stylize('[Function' + name + ']', 'special');
    }
    if (isRegExp(value)) {
      return ctx.stylize(RegExp.prototype.toString.call(value), 'regexp');
    }
    if (isDate(value)) {
      return ctx.stylize(Date.prototype.toString.call(value), 'date');
    }
    if (isError(value)) {
      return formatError(value);
    }
  }

  var base = '', array = false, braces = ['{', '}'];

  // Make Array say that they are Array
  if (isArray(value)) {
    array = true;
    braces = ['[', ']'];
  }

  // Make functions say that they are functions
  if (isFunction(value)) {
    var n = value.name ? ': ' + value.name : '';
    base = ' [Function' + n + ']';
  }

  // Make RegExps say that they are RegExps
  if (isRegExp(value)) {
    base = ' ' + RegExp.prototype.toString.call(value);
  }

  // Make dates with properties first say the date
  if (isDate(value)) {
    base = ' ' + Date.prototype.toUTCString.call(value);
  }

  // Make error with message first say the error
  if (isError(value)) {
    base = ' ' + formatError(value);
  }

  if (keys.length === 0 && (!array || value.length == 0)) {
    return braces[0] + base + braces[1];
  }

  if (recurseTimes < 0) {
    if (isRegExp(value)) {
      return ctx.stylize(RegExp.prototype.toString.call(value), 'regexp');
    } else {
      return ctx.stylize('[Object]', 'special');
    }
  }

  ctx.seen.push(value);

  var output;
  if (array) {
    output = formatArray(ctx, value, recurseTimes, visibleKeys, keys);
  } else {
    output = keys.map(function(key) {
      return formatProperty(ctx, value, recurseTimes, visibleKeys, key, array);
    });
  }

  ctx.seen.pop();

  return reduceToSingleString(output, base, braces);
}


function formatPrimitive(ctx, value) {
  if (isUndefined(value))
    return ctx.stylize('undefined', 'undefined');
  if (isString(value)) {
    var simple = '\'' + JSON.stringify(value).replace(/^"|"$/g, '')
                                             .replace(/'/g, "\\'")
                                             .replace(/\\"/g, '"') + '\'';
    return ctx.stylize(simple, 'string');
  }
  if (isNumber(value))
    return ctx.stylize('' + value, 'number');
  if (isBoolean(value))
    return ctx.stylize('' + value, 'boolean');
  // For some reason typeof null is "object", so special case here.
  if (isNull(value))
    return ctx.stylize('null', 'null');
}


function formatError(value) {
  return '[' + Error.prototype.toString.call(value) + ']';
}


function formatArray(ctx, value, recurseTimes, visibleKeys, keys) {
  var output = [];
  for (var i = 0, l = value.length; i < l; ++i) {
    if (hasOwnProperty(value, String(i))) {
      output.push(formatProperty(ctx, value, recurseTimes, visibleKeys,
          String(i), true));
    } else {
      output.push('');
    }
  }
  keys.forEach(function(key) {
    if (!key.match(/^\d+$/)) {
      output.push(formatProperty(ctx, value, recurseTimes, visibleKeys,
          key, true));
    }
  });
  return output;
}


function formatProperty(ctx, value, recurseTimes, visibleKeys, key, array) {
  var name, str, desc;
  desc = Object.getOwnPropertyDescriptor(value, key) || { value: value[key] };
  if (desc.get) {
    if (desc.set) {
      str = ctx.stylize('[Getter/Setter]', 'special');
    } else {
      str = ctx.stylize('[Getter]', 'special');
    }
  } else {
    if (desc.set) {
      str = ctx.stylize('[Setter]', 'special');
    }
  }
  if (!hasOwnProperty(visibleKeys, key)) {
    name = '[' + key + ']';
  }
  if (!str) {
    if (ctx.seen.indexOf(desc.value) < 0) {
      if (isNull(recurseTimes)) {
        str = formatValue(ctx, desc.value, null);
      } else {
        str = formatValue(ctx, desc.value, recurseTimes - 1);
      }
      if (str.indexOf('\n') > -1) {
        if (array) {
          str = str.split('\n').map(function(line) {
            return '  ' + line;
          }).join('\n').substr(2);
        } else {
          str = '\n' + str.split('\n').map(function(line) {
            return '   ' + line;
          }).join('\n');
        }
      }
    } else {
      str = ctx.stylize('[Circular]', 'special');
    }
  }
  if (isUndefined(name)) {
    if (array && key.match(/^\d+$/)) {
      return str;
    }
    name = JSON.stringify('' + key);
    if (name.match(/^"([a-zA-Z_][a-zA-Z_0-9]*)"$/)) {
      name = name.substr(1, name.length - 2);
      name = ctx.stylize(name, 'name');
    } else {
      name = name.replace(/'/g, "\\'")
                 .replace(/\\"/g, '"')
                 .replace(/(^"|"$)/g, "'");
      name = ctx.stylize(name, 'string');
    }
  }

  return name + ': ' + str;
}


function reduceToSingleString(output, base, braces) {
  var numLinesEst = 0;
  var length = output.reduce(function(prev, cur) {
    numLinesEst++;
    if (cur.indexOf('\n') >= 0) numLinesEst++;
    return prev + cur.replace(/\u001b\[\d\d?m/g, '').length + 1;
  }, 0);

  if (length > 60) {
    return braces[0] +
           (base === '' ? '' : base + '\n ') +
           ' ' +
           output.join(',\n  ') +
           ' ' +
           braces[1];
  }

  return braces[0] + base + ' ' + output.join(', ') + ' ' + braces[1];
}


// NOTE: These type checking functions intentionally don't use `instanceof`
// because it is fragile and can be easily faked with `Object.create()`.
function isArray(ar) {
  return Array.isArray(ar);
}
exports.isArray = isArray;

function isBoolean(arg) {
  return typeof arg === 'boolean';
}
exports.isBoolean = isBoolean;

function isNull(arg) {
  return arg === null;
}
exports.isNull = isNull;

function isNullOrUndefined(arg) {
  return arg == null;
}
exports.isNullOrUndefined = isNullOrUndefined;

function isNumber(arg) {
  return typeof arg === 'number';
}
exports.isNumber = isNumber;

function isString(arg) {
  return typeof arg === 'string';
}
exports.isString = isString;

function isSymbol(arg) {
  return typeof arg === 'symbol';
}
exports.isSymbol = isSymbol;

function isUndefined(arg) {
  return arg === void 0;
}
exports.isUndefined = isUndefined;

function isRegExp(re) {
  return isObject(re) && objectToString(re) === '[object RegExp]';
}
exports.isRegExp = isRegExp;

function isObject(arg) {
  return typeof arg === 'object' && arg !== null;
}
exports.isObject = isObject;

function isDate(d) {
  return isObject(d) && objectToString(d) === '[object Date]';
}
exports.isDate = isDate;

function isError(e) {
  return isObject(e) &&
      (objectToString(e) === '[object Error]' || e instanceof Error);
}
exports.isError = isError;

function isFunction(arg) {
  return typeof arg === 'function';
}
exports.isFunction = isFunction;

function isPrimitive(arg) {
  return arg === null ||
         typeof arg === 'boolean' ||
         typeof arg === 'number' ||
         typeof arg === 'string' ||
         typeof arg === 'symbol' ||  // ES6 symbol
         typeof arg === 'undefined';
}
exports.isPrimitive = isPrimitive;

exports.isBuffer = require('./support/isBuffer');

function objectToString(o) {
  return Object.prototype.toString.call(o);
}


function pad(n) {
  return n < 10 ? '0' + n.toString(10) : n.toString(10);
}


var months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep',
              'Oct', 'Nov', 'Dec'];

// 26 Feb 16:19:34
function timestamp() {
  var d = new Date();
  var time = [pad(d.getHours()),
              pad(d.getMinutes()),
              pad(d.getSeconds())].join(':');
  return [d.getDate(), months[d.getMonth()], time].join(' ');
}


// log is just a thin wrapper to console.log that prepends a timestamp
exports.log = function() {
  console.log('%s - %s', timestamp(), exports.format.apply(exports, arguments));
};


/**
 * Inherit the prototype methods from one constructor into another.
 *
 * The Function.prototype.inherits from lang.js rewritten as a standalone
 * function (not on Function.prototype). NOTE: If this file is to be loaded
 * during bootstrapping this function needs to be rewritten using some native
 * functions as prototype setup using normal JavaScript does not work as
 * expected during bootstrapping (see mirror.js in r114903).
 *
 * @param {function} ctor Constructor function which needs to inherit the
 *     prototype.
 * @param {function} superCtor Constructor function to inherit prototype from.
 */
exports.inherits = require('inherits');

exports._extend = function(origin, add) {
  // Don't do anything if add isn't an object
  if (!add || !isObject(add)) return origin;

  var keys = Object.keys(add);
  var i = keys.length;
  while (i--) {
    origin[keys[i]] = add[keys[i]];
  }
  return origin;
};

function hasOwnProperty(obj, prop) {
  return Object.prototype.hasOwnProperty.call(obj, prop);
}

}).call(this,require('_process'),typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./support/isBuffer":"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/util/support/isBufferBrowser.js","_process":"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/process/browser.js","inherits":"/usr/local/lib/node_modules/watchify/node_modules/browserify/node_modules/inherits/inherits_browser.js"}]},{},["/Users/vallette/Desktop/colony/src/index.js"])
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi91c3IvbG9jYWwvbGliL25vZGVfbW9kdWxlcy93YXRjaGlmeS9ub2RlX21vZHVsZXMvYnJvd3NlcmlmeS9ub2RlX21vZHVsZXMvYnJvd3Nlci1wYWNrL19wcmVsdWRlLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvaWNoLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXNjYWxlL25vZGVfbW9kdWxlcy90d28tc3VtL3R3by1zdW0uanMiLCIvVXNlcnMvdmFsbGV0dGUvRGVza3RvcC9jb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL25vZGVfbW9kdWxlcy9yb2J1c3Qtc2NhbGUvcm9idXN0LXNjYWxlLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1YnRyYWN0L3JvYnVzdC1kaWZmLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvcm9idXN0LXN1bS9yb2J1c3Qtc3VtLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3JvYnVzdC1vcmllbnRhdGlvbi9ub2RlX21vZHVsZXMvdHdvLXByb2R1Y3QvdHdvLXByb2R1Y3QuanMiLCIvVXNlcnMvdmFsbGV0dGUvRGVza3RvcC9jb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvcm9idXN0LW9yaWVudGF0aW9uL29yaWVudGF0aW9uLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvaW5jcmVtZW50YWwtY29udmV4LWh1bGwvbm9kZV9tb2R1bGVzL3NpbXBsaWNpYWwtY29tcGxleC9ub2RlX21vZHVsZXMvYml0LXR3aWRkbGUvdHdpZGRsZS5qcyIsIi9Vc2Vycy92YWxsZXR0ZS9EZXNrdG9wL2NvbG9ueS9ub2RlX21vZHVsZXMvZGVsYXVuYXktdHJpYW5ndWxhdGUvbm9kZV9tb2R1bGVzL2luY3JlbWVudGFsLWNvbnZleC1odWxsL25vZGVfbW9kdWxlcy9zaW1wbGljaWFsLWNvbXBsZXgvbm9kZV9tb2R1bGVzL3VuaW9uLWZpbmQvaW5kZXguanMiLCIvVXNlcnMvdmFsbGV0dGUvRGVza3RvcC9jb2xvbnkvbm9kZV9tb2R1bGVzL2RlbGF1bmF5LXRyaWFuZ3VsYXRlL25vZGVfbW9kdWxlcy9pbmNyZW1lbnRhbC1jb252ZXgtaHVsbC9ub2RlX21vZHVsZXMvc2ltcGxpY2lhbC1jb21wbGV4L3RvcG9sb2d5LmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS9ub2RlX21vZHVsZXMvdW5pcS91bmlxLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9kZWxhdW5heS10cmlhbmd1bGF0ZS90cmlhbmd1bGF0ZS5qcyIsIi9Vc2Vycy92YWxsZXR0ZS9EZXNrdG9wL2NvbG9ueS9ub2RlX21vZHVsZXMvZ2FtZS1zaGVsbC9saWIvaHJ0aW1lLXBvbHlmaWxsLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9nYW1lLXNoZWxsL2xpYi9tb3VzZXdoZWVsLXBvbHlmaWxsLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9nYW1lLXNoZWxsL2xpYi9yYWYtcG9seWZpbGwuanMiLCIvVXNlcnMvdmFsbGV0dGUvRGVza3RvcC9jb2xvbnkvbm9kZV9tb2R1bGVzL2dhbWUtc2hlbGwvbm9kZV9tb2R1bGVzL2JpbmFyeS1zZWFyY2gtYm91bmRzL3NlYXJjaC1ib3VuZHMuanMiLCIvVXNlcnMvdmFsbGV0dGUvRGVza3RvcC9jb2xvbnkvbm9kZV9tb2R1bGVzL2dhbWUtc2hlbGwvbm9kZV9tb2R1bGVzL2RvbXJlYWR5L3JlYWR5LmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9nYW1lLXNoZWxsL25vZGVfbW9kdWxlcy9pbnZlcnQtaGFzaC9pbnZlcnQuanMiLCIvVXNlcnMvdmFsbGV0dGUvRGVza3RvcC9jb2xvbnkvbm9kZV9tb2R1bGVzL2dhbWUtc2hlbGwvbm9kZV9tb2R1bGVzL2lvdGEtYXJyYXkvaW90YS5qcyIsIi9Vc2Vycy92YWxsZXR0ZS9EZXNrdG9wL2NvbG9ueS9ub2RlX21vZHVsZXMvZ2FtZS1zaGVsbC9ub2RlX21vZHVsZXMvdW5pcS91bmlxLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9nYW1lLXNoZWxsL25vZGVfbW9kdWxlcy92a2V5L2luZGV4LmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9nYW1lLXNoZWxsL3NoZWxsLmpzIiwiL1VzZXJzL3ZhbGxldHRlL0Rlc2t0b3AvY29sb255L25vZGVfbW9kdWxlcy9oYXJtb255LWNvbGxlY3Rpb25zL2hhcm1vbnktY29sbGVjdGlvbnMuanMiLCIvVXNlcnMvdmFsbGV0dGUvRGVza3RvcC9jb2xvbnkvc3JjL2luZGV4LmpzIiwiL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL3dhdGNoaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9ldmVudHMvZXZlbnRzLmpzIiwiL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL3dhdGNoaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9pbmhlcml0cy9pbmhlcml0c19icm93c2VyLmpzIiwiL3Vzci9sb2NhbC9saWIvbm9kZV9tb2R1bGVzL3dhdGNoaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9wcm9jZXNzL2Jyb3dzZXIuanMiLCIvdXNyL2xvY2FsL2xpYi9ub2RlX21vZHVsZXMvd2F0Y2hpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL3V0aWwvc3VwcG9ydC9pc0J1ZmZlckJyb3dzZXIuanMiLCIvdXNyL2xvY2FsL2xpYi9ub2RlX21vZHVsZXMvd2F0Y2hpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL3V0aWwvdXRpbC5qcyJdLCJuYW1lcyI6W10sIm1hcHBpbmdzIjoiQUFBQTtBQ0FBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDN2JBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDaEJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDakRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMzSkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzNKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDaENBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzdMQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM1TUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDekRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3RWQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN6REE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzlKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDWEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMxREE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM1QkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNURBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNUJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ1pBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDVkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDekRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDeElBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOXNCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMzeEJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOU9BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQy9TQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdkJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQy9EQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDTEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EiLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbiBlKHQsbixyKXtmdW5jdGlvbiBzKG8sdSl7aWYoIW5bb10pe2lmKCF0W29dKXt2YXIgYT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2lmKCF1JiZhKXJldHVybiBhKG8sITApO2lmKGkpcmV0dXJuIGkobywhMCk7dmFyIGY9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitvK1wiJ1wiKTt0aHJvdyBmLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsZn12YXIgbD1uW29dPXtleHBvcnRzOnt9fTt0W29dWzBdLmNhbGwobC5leHBvcnRzLGZ1bmN0aW9uKGUpe3ZhciBuPXRbb11bMV1bZV07cmV0dXJuIHMobj9uOmUpfSxsLGwuZXhwb3J0cyxlLHQsbixyKX1yZXR1cm4gbltvXS5leHBvcnRzfXZhciBpPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7Zm9yKHZhciBvPTA7bzxyLmxlbmd0aDtvKyspcyhyW29dKTtyZXR1cm4gc30pIiwiXCJ1c2Ugc3RyaWN0XCJcblxuLy9IaWdoIGxldmVsIGlkZWE6XG4vLyAxLiBVc2UgQ2xhcmtzb24ncyBpbmNyZW1lbnRhbCBjb25zdHJ1Y3Rpb24gdG8gZmluZCBjb252ZXggaHVsbFxuLy8gMi4gUG9pbnQgbG9jYXRpb24gaW4gdHJpYW5ndWxhdGlvbiBieSBqdW1wIGFuZCB3YWxrXG5cbm1vZHVsZS5leHBvcnRzID0gaW5jcmVtZW50YWxDb252ZXhIdWxsXG5cbnZhciBvcmllbnQgPSByZXF1aXJlKFwicm9idXN0LW9yaWVudGF0aW9uXCIpXG52YXIgY29tcGFyZUNlbGwgPSByZXF1aXJlKFwic2ltcGxpY2lhbC1jb21wbGV4XCIpLmNvbXBhcmVDZWxsc1xuXG5mdW5jdGlvbiBjb21wYXJlSW50KGEsIGIpIHtcbiAgcmV0dXJuIGEgLSBiXG59XG5cbmZ1bmN0aW9uIFNpbXBsZXgodmVydGljZXMsIGFkamFjZW50LCBib3VuZGFyeSkge1xuICB0aGlzLnZlcnRpY2VzID0gdmVydGljZXNcbiAgdGhpcy5hZGphY2VudCA9IGFkamFjZW50XG4gIHRoaXMuYm91bmRhcnkgPSBib3VuZGFyeVxuICB0aGlzLmxhc3RWaXNpdGVkID0gLTFcbn1cblxuU2ltcGxleC5wcm90b3R5cGUuZmxpcCA9IGZ1bmN0aW9uKCkge1xuICB2YXIgdCA9IHRoaXMudmVydGljZXNbMF1cbiAgdGhpcy52ZXJ0aWNlc1swXSA9IHRoaXMudmVydGljZXNbMV1cbiAgdGhpcy52ZXJ0aWNlc1sxXSA9IHRcbiAgdmFyIHUgPSB0aGlzLmFkamFjZW50WzBdXG4gIHRoaXMuYWRqYWNlbnRbMF0gPSB0aGlzLmFkamFjZW50WzFdXG4gIHRoaXMuYWRqYWNlbnRbMV0gPSB1XG59XG5cbmZ1bmN0aW9uIEdsdWVGYWNldCh2ZXJ0aWNlcywgY2VsbCwgaW5kZXgpIHtcbiAgdGhpcy52ZXJ0aWNlcyA9IHZlcnRpY2VzXG4gIHRoaXMuY2VsbCA9IGNlbGxcbiAgdGhpcy5pbmRleCA9IGluZGV4XG59XG5cbmZ1bmN0aW9uIGNvbXBhcmVHbHVlKGEsIGIpIHtcbiAgcmV0dXJuIGNvbXBhcmVDZWxsKGEudmVydGljZXMsIGIudmVydGljZXMpXG59XG5cbmZ1bmN0aW9uIGJha2VPcmllbnQoZCkge1xuICB2YXIgY29kZSA9IFtcImZ1bmN0aW9uIG9yaWVudCgpe3ZhciB0dXBsZT10aGlzLnR1cGxlO3JldHVybiB0ZXN0KFwiXVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgaWYoaSA+IDApIHtcbiAgICAgIGNvZGUucHVzaChcIixcIilcbiAgICB9XG4gICAgY29kZS5wdXNoKFwidHVwbGVbXCIsIGksIFwiXVwiKVxuICB9XG4gIGNvZGUucHVzaChcIil9cmV0dXJuIG9yaWVudFwiKVxuICB2YXIgcHJvYyA9IG5ldyBGdW5jdGlvbihcInRlc3RcIiwgY29kZS5qb2luKFwiXCIpKVxuICB2YXIgdGVzdCA9IG9yaWVudFtkKzFdXG4gIGlmKCF0ZXN0KSB7XG4gICAgdGVzdCA9IG9yaWVudFxuICB9XG4gIHJldHVybiBwcm9jKHRlc3QpXG59XG5cbnZhciBCQUtFRCA9IFtdXG5cbmZ1bmN0aW9uIFRyaWFuZ3VsYXRpb24oZGltZW5zaW9uLCB2ZXJ0aWNlcywgc2ltcGxpY2VzKSB7XG4gIHRoaXMuZGltZW5zaW9uID0gZGltZW5zaW9uXG4gIHRoaXMudmVydGljZXMgPSB2ZXJ0aWNlc1xuICB0aGlzLnNpbXBsaWNlcyA9IHNpbXBsaWNlc1xuICB0aGlzLmludGVyaW9yID0gc2ltcGxpY2VzLmZpbHRlcihmdW5jdGlvbihjKSB7XG4gICAgcmV0dXJuICFjLmJvdW5kYXJ5XG4gIH0pXG5cbiAgdGhpcy50dXBsZSA9IG5ldyBBcnJheShkaW1lbnNpb24rMSlcbiAgZm9yKHZhciBpPTA7IGk8PWRpbWVuc2lvbjsgKytpKSB7XG4gICAgdGhpcy50dXBsZVtpXSA9IHRoaXMudmVydGljZXNbaV1cbiAgfVxuXG4gIHZhciBvID0gQkFLRURbZGltZW5zaW9uXVxuICBpZighbykge1xuICAgIG8gPSBCQUtFRFtkaW1lbnNpb25dID0gYmFrZU9yaWVudChkaW1lbnNpb24pXG4gIH1cbiAgdGhpcy5vcmllbnQgPSBvXG59XG5cbnZhciBwcm90byA9IFRyaWFuZ3VsYXRpb24ucHJvdG90eXBlXG5cbi8vRGVnZW5lcmF0ZSBzaXR1YXRpb24gd2hlcmUgd2UgYXJlIG9uIGJvdW5kYXJ5LCBidXQgY29wbGFuYXIgdG8gZmFjZVxucHJvdG8uaGFuZGxlQm91bmRhcnlEZWdlbmVyYWN5ID0gZnVuY3Rpb24oY2VsbCwgcG9pbnQpIHtcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgbiA9IHRoaXMudmVydGljZXMubGVuZ3RoIC0gMVxuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG4gIHZhciB2ZXJ0cyA9IHRoaXMudmVydGljZXNcblxuICAvL0R1bWIgc29sdXRpb246IEp1c3QgZG8gZGZzIGZyb20gYm91bmRhcnkgY2VsbCB1bnRpbCB3ZSBmaW5kIGFueSBwZWFrLCBvciB0ZXJtaW5hdGVcbiAgdmFyIHRvVmlzaXQgPSBbIGNlbGwgXVxuICBjZWxsLmxhc3RWaXNpdGVkID0gLW5cbiAgd2hpbGUodG9WaXNpdC5sZW5ndGggPiAwKSB7XG4gICAgY2VsbCA9IHRvVmlzaXQucG9wKClcbiAgICB2YXIgY2VsbFZlcnRzID0gY2VsbC52ZXJ0aWNlc1xuICAgIHZhciBjZWxsQWRqID0gY2VsbC5hZGphY2VudFxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIHZhciBuZWlnaGJvciA9IGNlbGxBZGpbaV1cbiAgICAgIGlmKCFuZWlnaGJvci5ib3VuZGFyeSB8fCBuZWlnaGJvci5sYXN0VmlzaXRlZCA8PSAtbikge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgdmFyIG52ID0gbmVpZ2hib3IudmVydGljZXNcbiAgICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgICAgdmFyIHZ2ID0gbnZbal1cbiAgICAgICAgaWYodnYgPCAwKSB7XG4gICAgICAgICAgdHVwbGVbal0gPSBwb2ludFxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgIHR1cGxlW2pdID0gdmVydHNbdnZdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHZhciBvID0gdGhpcy5vcmllbnQoKVxuICAgICAgaWYobyA+IDApIHtcbiAgICAgICAgcmV0dXJuIG5laWdoYm9yXG4gICAgICB9XG4gICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IC1uXG4gICAgICBpZihvID09PSAwKSB7XG4gICAgICAgIHRvVmlzaXQucHVzaChuZWlnaGJvcilcbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgcmV0dXJuIG51bGxcbn1cblxucHJvdG8ud2FsayA9IGZ1bmN0aW9uKHBvaW50LCByYW5kb20pIHtcbiAgLy9BbGlhcyBsb2NhbCBwcm9wZXJ0aWVzXG4gIHZhciBuID0gdGhpcy52ZXJ0aWNlcy5sZW5ndGggLSAxXG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG5cbiAgLy9Db21wdXRlIGluaXRpYWwganVtcCBjZWxsXG4gIHZhciBpbml0SW5kZXggPSByYW5kb20gPyAodGhpcy5pbnRlcmlvci5sZW5ndGggKiBNYXRoLnJhbmRvbSgpKXwwIDogKHRoaXMuaW50ZXJpb3IubGVuZ3RoLTEpXG4gIHZhciBjZWxsID0gdGhpcy5pbnRlcmlvclsgaW5pdEluZGV4IF1cblxuICAvL1N0YXJ0IHdhbGtpbmdcbm91dGVyTG9vcDpcbiAgd2hpbGUoIWNlbGwuYm91bmRhcnkpIHtcbiAgICB2YXIgY2VsbFZlcnRzID0gY2VsbC52ZXJ0aWNlc1xuICAgIHZhciBjZWxsQWRqID0gY2VsbC5hZGphY2VudFxuXG4gICAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgICAgdHVwbGVbaV0gPSB2ZXJ0c1tjZWxsVmVydHNbaV1dXG4gICAgfVxuICAgIGNlbGwubGFzdFZpc2l0ZWQgPSBuXG5cbiAgICAvL0ZpbmQgZmFydGhlc3QgYWRqYWNlbnQgY2VsbFxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIHZhciBuZWlnaGJvciA9IGNlbGxBZGpbaV1cbiAgICAgIGlmKG5laWdoYm9yLmxhc3RWaXNpdGVkID49IG4pIHtcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIHZhciBwcmV2ID0gdHVwbGVbaV1cbiAgICAgIHR1cGxlW2ldID0gcG9pbnRcbiAgICAgIHZhciBvID0gdGhpcy5vcmllbnQoKVxuICAgICAgdHVwbGVbaV0gPSBwcmV2XG4gICAgICBpZihvIDwgMCkge1xuICAgICAgICBjZWxsID0gbmVpZ2hib3JcbiAgICAgICAgY29udGludWUgb3V0ZXJMb29wXG4gICAgICB9IGVsc2Uge1xuICAgICAgICBpZighbmVpZ2hib3IuYm91bmRhcnkpIHtcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IG5cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICBuZWlnaGJvci5sYXN0VmlzaXRlZCA9IC1uXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gICAgcmV0dXJuXG4gIH1cblxuICByZXR1cm4gY2VsbFxufVxuXG5wcm90by5hZGRQZWFrcyA9IGZ1bmN0aW9uKHBvaW50LCBjZWxsKSB7XG4gIHZhciBuID0gdGhpcy52ZXJ0aWNlcy5sZW5ndGggLSAxXG4gIHZhciBkID0gdGhpcy5kaW1lbnNpb25cbiAgdmFyIHZlcnRzID0gdGhpcy52ZXJ0aWNlc1xuICB2YXIgdHVwbGUgPSB0aGlzLnR1cGxlXG4gIHZhciBpbnRlcmlvciA9IHRoaXMuaW50ZXJpb3JcbiAgdmFyIHNpbXBsaWNlcyA9IHRoaXMuc2ltcGxpY2VzXG5cbiAgLy9XYWxraW5nIGZpbmlzaGVkIGF0IGJvdW5kYXJ5LCB0aW1lIHRvIGFkZCBwZWFrc1xuICB2YXIgdG92aXNpdCA9IFsgY2VsbCBdXG5cbiAgLy9TdHJldGNoIGluaXRpYWwgYm91bmRhcnkgY2VsbCBpbnRvIGEgcGVha1xuICBjZWxsLmxhc3RWaXNpdGVkID0gblxuICBjZWxsLnZlcnRpY2VzW2NlbGwudmVydGljZXMuaW5kZXhPZigtMSldID0gblxuICBjZWxsLmJvdW5kYXJ5ID0gZmFsc2VcbiAgaW50ZXJpb3IucHVzaChjZWxsKVxuXG4gIC8vUmVjb3JkIGEgbGlzdCBvZiBhbGwgbmV3IGJvdW5kYXJpZXMgY3JlYXRlZCBieSBhZGRlZCBwZWFrcyBzbyB3ZSBjYW4gZ2x1ZSB0aGVtIHRvZ2V0aGVyIHdoZW4gd2UgYXJlIGFsbCBkb25lXG4gIHZhciBnbHVlRmFjZXRzID0gW11cblxuICAvL0RvIGEgdHJhdmVyc2FsIG9mIHRoZSBib3VuZGFyeSB3YWxraW5nIG91dHdhcmQgZnJvbSBzdGFydGluZyBwZWFrXG4gIHdoaWxlKHRvdmlzaXQubGVuZ3RoID4gMCkge1xuICAgIC8vUG9wIG9mZiBwZWFrIGFuZCB3YWxrIG92ZXIgYWRqYWNlbnQgY2VsbHNcbiAgICB2YXIgY2VsbCA9IHRvdmlzaXQucG9wKClcbiAgICB2YXIgY2VsbFZlcnRzID0gY2VsbC52ZXJ0aWNlc1xuICAgIHZhciBjZWxsQWRqID0gY2VsbC5hZGphY2VudFxuICAgIHZhciBpbmRleE9mTiA9IGNlbGxWZXJ0cy5pbmRleE9mKG4pXG4gICAgaWYoaW5kZXhPZk4gPCAwKSB7XG4gICAgICBjb250aW51ZVxuICAgIH1cblxuICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgIGlmKGkgPT09IGluZGV4T2ZOKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG5cbiAgICAgIC8vRm9yIGVhY2ggYm91bmRhcnkgbmVpZ2hib3Igb2YgdGhlIGNlbGxcbiAgICAgIHZhciBuZWlnaGJvciA9IGNlbGxBZGpbaV1cbiAgICAgIGlmKCFuZWlnaGJvci5ib3VuZGFyeSB8fCBuZWlnaGJvci5sYXN0VmlzaXRlZCA+PSBuKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG5cbiAgICAgIHZhciBudiA9IG5laWdoYm9yLnZlcnRpY2VzXG5cbiAgICAgIC8vVGVzdCBpZiBuZWlnaGJvciBpcyBhIHBlYWtcbiAgICAgIGlmKG5laWdoYm9yLmxhc3RWaXNpdGVkICE9PSAtbikgeyAgICAgIFxuICAgICAgICAvL0NvbXB1dGUgb3JpZW50YXRpb24gb2YgcCByZWxhdGl2ZSB0byBlYWNoIGJvdW5kYXJ5IHBlYWtcbiAgICAgICAgdmFyIGluZGV4T2ZOZWcxID0gMFxuICAgICAgICBmb3IodmFyIGo9MDsgajw9ZDsgKytqKSB7XG4gICAgICAgICAgaWYobnZbal0gPCAwKSB7XG4gICAgICAgICAgICBpbmRleE9mTmVnMSA9IGpcbiAgICAgICAgICAgIHR1cGxlW2pdID0gcG9pbnRcbiAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdHVwbGVbal0gPSB2ZXJ0c1tudltqXV1cbiAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgdmFyIG8gPSB0aGlzLm9yaWVudCgpXG5cbiAgICAgICAgLy9UZXN0IGlmIG5laWdoYm9yIGNlbGwgaXMgYWxzbyBhIHBlYWtcbiAgICAgICAgaWYobyA+IDApIHtcbiAgICAgICAgICBudltpbmRleE9mTmVnMV0gPSBuXG4gICAgICAgICAgbmVpZ2hib3IuYm91bmRhcnkgPSBmYWxzZVxuICAgICAgICAgIGludGVyaW9yLnB1c2gobmVpZ2hib3IpXG4gICAgICAgICAgdG92aXNpdC5wdXNoKG5laWdoYm9yKVxuICAgICAgICAgIG5laWdoYm9yLmxhc3RWaXNpdGVkID0gblxuICAgICAgICAgIGNvbnRpbnVlXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgbmVpZ2hib3IubGFzdFZpc2l0ZWQgPSAtblxuICAgICAgICB9XG4gICAgICB9XG5cbiAgICAgIHZhciBuYSA9IG5laWdoYm9yLmFkamFjZW50XG5cbiAgICAgIC8vT3RoZXJ3aXNlLCByZXBsYWNlIG5laWdoYm9yIHdpdGggbmV3IGZhY2VcbiAgICAgIHZhciB2dmVydHMgPSBjZWxsVmVydHMuc2xpY2UoKVxuICAgICAgdmFyIHZhZGogPSBjZWxsQWRqLnNsaWNlKClcbiAgICAgIHZhciBuY2VsbCA9IG5ldyBTaW1wbGV4KHZ2ZXJ0cywgdmFkaiwgdHJ1ZSlcbiAgICAgIHNpbXBsaWNlcy5wdXNoKG5jZWxsKVxuXG4gICAgICAvL0Nvbm5lY3QgdG8gbmVpZ2hib3JcbiAgICAgIHZhciBvcHBvc2l0ZSA9IG5hLmluZGV4T2YoY2VsbClcbiAgICAgIGlmKG9wcG9zaXRlIDwgMCkge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbmFbb3Bwb3NpdGVdID0gbmNlbGxcbiAgICAgIHZhZGpbaW5kZXhPZk5dID0gbmVpZ2hib3JcblxuICAgICAgLy9Db25uZWN0IHRvIGNlbGxcbiAgICAgIHZ2ZXJ0c1tpXSA9IC0xXG4gICAgICB2YWRqW2ldID0gY2VsbFxuICAgICAgY2VsbEFkaltpXSA9IG5jZWxsXG5cbiAgICAgIC8vRmxpcCBmYWNldFxuICAgICAgbmNlbGwuZmxpcCgpXG5cbiAgICAgIC8vQWRkIHRvIGdsdWUgbGlzdFxuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICB2YXIgdXUgPSB2dmVydHNbal1cbiAgICAgICAgaWYodXUgPCAwIHx8IHV1ID09PSBuKSB7XG4gICAgICAgICAgY29udGludWVcbiAgICAgICAgfVxuICAgICAgICB2YXIgbmZhY2UgPSBuZXcgQXJyYXkoZC0xKVxuICAgICAgICB2YXIgbnB0ciA9IDBcbiAgICAgICAgZm9yKHZhciBrPTA7IGs8PWQ7ICsraykge1xuICAgICAgICAgIHZhciB2diA9IHZ2ZXJ0c1trXVxuICAgICAgICAgIGlmKHZ2IDwgMCB8fCBrID09PSBqKSB7XG4gICAgICAgICAgICBjb250aW51ZVxuICAgICAgICAgIH1cbiAgICAgICAgICBuZmFjZVtucHRyKytdID0gdnZcbiAgICAgICAgfVxuICAgICAgICBnbHVlRmFjZXRzLnB1c2gobmV3IEdsdWVGYWNldChuZmFjZSwgbmNlbGwsIGopKVxuICAgICAgfVxuICAgIH1cbiAgfVxuXG4gIC8vR2x1ZSBib3VuZGFyeSBmYWNldHMgdG9nZXRoZXJcbiAgZ2x1ZUZhY2V0cy5zb3J0KGNvbXBhcmVHbHVlKVxuXG4gIGZvcih2YXIgaT0wOyBpKzE8Z2x1ZUZhY2V0cy5sZW5ndGg7IGkrPTIpIHtcbiAgICB2YXIgYSA9IGdsdWVGYWNldHNbaV1cbiAgICB2YXIgYiA9IGdsdWVGYWNldHNbaSsxXVxuICAgIHZhciBhaSA9IGEuaW5kZXhcbiAgICB2YXIgYmkgPSBiLmluZGV4XG4gICAgaWYoYWkgPCAwIHx8IGJpIDwgMCkge1xuICAgICAgY29udGludWVcbiAgICB9XG4gICAgYS5jZWxsLmFkamFjZW50W2EuaW5kZXhdID0gYi5jZWxsXG4gICAgYi5jZWxsLmFkamFjZW50W2IuaW5kZXhdID0gYS5jZWxsXG4gIH1cbn1cblxucHJvdG8uaW5zZXJ0ID0gZnVuY3Rpb24ocG9pbnQsIHJhbmRvbSkge1xuICAvL0FkZCBwb2ludFxuICB2YXIgdmVydHMgPSB0aGlzLnZlcnRpY2VzXG4gIHZlcnRzLnB1c2gocG9pbnQpXG5cbiAgdmFyIGNlbGwgPSB0aGlzLndhbGsocG9pbnQsIHJhbmRvbSlcbiAgaWYoIWNlbGwpIHtcbiAgICByZXR1cm5cbiAgfVxuXG4gIC8vQWxpYXMgbG9jYWwgcHJvcGVydGllc1xuICB2YXIgZCA9IHRoaXMuZGltZW5zaW9uXG4gIHZhciB0dXBsZSA9IHRoaXMudHVwbGVcblxuICAvL0RlZ2VuZXJhdGUgY2FzZTogSWYgcG9pbnQgaXMgY29wbGFuYXIgdG8gY2VsbCwgdGhlbiB3YWxrIHVudGlsIHdlIGZpbmQgYSBub24tZGVnZW5lcmF0ZSBib3VuZGFyeVxuICBmb3IodmFyIGk9MDsgaTw9ZDsgKytpKSB7XG4gICAgdmFyIHZ2ID0gY2VsbC52ZXJ0aWNlc1tpXVxuICAgIGlmKHZ2IDwgMCkge1xuICAgICAgdHVwbGVbaV0gPSBwb2ludFxuICAgIH0gZWxzZSB7XG4gICAgICB0dXBsZVtpXSA9IHZlcnRzW3Z2XVxuICAgIH1cbiAgfVxuICB2YXIgbyA9IHRoaXMub3JpZW50KHR1cGxlKVxuICBpZihvIDwgMCkge1xuICAgIHJldHVyblxuICB9IGVsc2UgaWYobyA9PT0gMCkge1xuICAgIGNlbGwgPSB0aGlzLmhhbmRsZUJvdW5kYXJ5RGVnZW5lcmFjeShjZWxsLCBwb2ludClcbiAgICBpZighY2VsbCkge1xuICAgICAgcmV0dXJuXG4gICAgfVxuICB9XG5cbiAgLy9BZGQgcGVha3NcbiAgdGhpcy5hZGRQZWFrcyhwb2ludCwgY2VsbClcbn1cblxuLy9FeHRyYWN0IGFsbCBib3VuZGFyeSBjZWxsc1xucHJvdG8uYm91bmRhcnkgPSBmdW5jdGlvbigpIHtcbiAgdmFyIGQgPSB0aGlzLmRpbWVuc2lvblxuICB2YXIgYm91bmRhcnkgPSBbXVxuICB2YXIgY2VsbHMgPSB0aGlzLnNpbXBsaWNlc1xuICB2YXIgbmMgPSBjZWxscy5sZW5ndGhcbiAgZm9yKHZhciBpPTA7IGk8bmM7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBpZihjLmJvdW5kYXJ5KSB7XG4gICAgICB2YXIgYmNlbGwgPSBuZXcgQXJyYXkoZClcbiAgICAgIHZhciBjdiA9IGMudmVydGljZXNcbiAgICAgIHZhciBwdHIgPSAwXG4gICAgICB2YXIgcGFyaXR5ID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICBpZihjdltqXSA+PSAwKSB7XG4gICAgICAgICAgYmNlbGxbcHRyKytdID0gY3Zbal1cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICBwYXJpdHkgPSBqJjFcbiAgICAgICAgfVxuICAgICAgfVxuICAgICAgaWYocGFyaXR5ID09PSAoZCYxKSkge1xuICAgICAgICB2YXIgdCA9IGJjZWxsWzBdXG4gICAgICAgIGJjZWxsWzBdID0gYmNlbGxbMV1cbiAgICAgICAgYmNlbGxbMV0gPSB0XG4gICAgICB9XG4gICAgICBib3VuZGFyeS5wdXNoKGJjZWxsKVxuICAgIH1cbiAgfVxuICByZXR1cm4gYm91bmRhcnlcbn1cblxuZnVuY3Rpb24gaW5jcmVtZW50YWxDb252ZXhIdWxsKHBvaW50cywgcmFuZG9tU2VhcmNoKSB7XG4gIHZhciBuID0gcG9pbnRzLmxlbmd0aFxuICBpZihuID09PSAwKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKFwiTXVzdCBoYXZlIGF0IGxlYXN0IGQrMSBwb2ludHNcIilcbiAgfVxuICB2YXIgZCA9IHBvaW50c1swXS5sZW5ndGhcbiAgaWYobiA8PSBkKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKFwiTXVzdCBpbnB1dCBhdCBsZWFzdCBkKzEgcG9pbnRzXCIpXG4gIH1cblxuICAvL0ZJWE1FOiBUaGlzIGNvdWxkIGJlIGRlZ2VuZXJhdGUsIGJ1dCBuZWVkIHRvIHNlbGVjdCBkKzEgbm9uLWNvcGxhbmFyIHBvaW50cyB0byBib290c3RyYXAgcHJvY2Vzc1xuICB2YXIgaW5pdGlhbFNpbXBsZXggPSBwb2ludHMuc2xpY2UoMCwgZCsxKVxuXG4gIC8vTWFrZSBzdXJlIGluaXRpYWwgc2ltcGxleCBpcyBwb3NpdGl2ZWx5IG9yaWVudGVkXG4gIHZhciBvID0gb3JpZW50LmFwcGx5KHZvaWQgMCwgaW5pdGlhbFNpbXBsZXgpXG4gIGlmKG8gPT09IDApIHtcbiAgICB0aHJvdyBuZXcgRXJyb3IoXCJJbnB1dCBub3QgaW4gZ2VuZXJhbCBwb3NpdGlvblwiKVxuICB9XG4gIHZhciBpbml0aWFsQ29vcmRzID0gbmV3IEFycmF5KGQrMSlcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIGluaXRpYWxDb29yZHNbaV0gPSBpXG4gIH1cbiAgaWYobyA8IDApIHtcbiAgICBpbml0aWFsQ29vcmRzWzBdID0gMVxuICAgIGluaXRpYWxDb29yZHNbMV0gPSAwXG4gIH1cblxuICAvL0NyZWF0ZSBpbml0aWFsIHRvcG9sb2dpY2FsIGluZGV4LCBnbHVlIHBvaW50ZXJzIHRvZ2V0aGVyIChraW5kIG9mIG1lc3N5KVxuICB2YXIgaW5pdGlhbENlbGwgPSBuZXcgU2ltcGxleChpbml0aWFsQ29vcmRzLCBuZXcgQXJyYXkoZCsxKSwgZmFsc2UpXG4gIHZhciBib3VuZGFyeSA9IGluaXRpYWxDZWxsLmFkamFjZW50XG4gIHZhciBsaXN0ID0gbmV3IEFycmF5KGQrMilcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB2ZXJ0cyA9IGluaXRpYWxDb29yZHMuc2xpY2UoKVxuICAgIGZvcih2YXIgaj0wOyBqPD1kOyArK2opIHtcbiAgICAgIGlmKGogPT09IGkpIHtcbiAgICAgICAgdmVydHNbal0gPSAtMVxuICAgICAgfVxuICAgIH1cbiAgICB2YXIgdCA9IHZlcnRzWzBdXG4gICAgdmVydHNbMF0gPSB2ZXJ0c1sxXVxuICAgIHZlcnRzWzFdID0gdFxuICAgIHZhciBjZWxsID0gbmV3IFNpbXBsZXgodmVydHMsIG5ldyBBcnJheShkKzEpLCB0cnVlKVxuICAgIGJvdW5kYXJ5W2ldID0gY2VsbFxuICAgIGxpc3RbaV0gPSBjZWxsXG4gIH1cbiAgbGlzdFtkKzFdID0gaW5pdGlhbENlbGxcbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB2ZXJ0cyA9IGJvdW5kYXJ5W2ldLnZlcnRpY2VzXG4gICAgdmFyIGFkaiA9IGJvdW5kYXJ5W2ldLmFkamFjZW50XG4gICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgdmFyIHYgPSB2ZXJ0c1tqXVxuICAgICAgaWYodiA8IDApIHtcbiAgICAgICAgYWRqW2pdID0gaW5pdGlhbENlbGxcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGZvcih2YXIgaz0wOyBrPD1kOyArK2spIHtcbiAgICAgICAgaWYoYm91bmRhcnlba10udmVydGljZXMuaW5kZXhPZih2KSA8IDApIHtcbiAgICAgICAgICBhZGpbal0gPSBib3VuZGFyeVtrXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgfVxuICB9XG5cbiAgLy9Jbml0aWFsaXplIHRyaWFuZ2xlc1xuICB2YXIgdHJpYW5nbGVzID0gbmV3IFRyaWFuZ3VsYXRpb24oZCwgaW5pdGlhbFNpbXBsZXgsIGxpc3QpXG5cbiAgLy9JbnNlcnQgcmVtYWluaW5nIHBvaW50c1xuICB2YXIgdXNlUmFuZG9tID0gISFyYW5kb21TZWFyY2hcbiAgZm9yKHZhciBpPWQrMTsgaTxuOyArK2kpIHtcbiAgICB0cmlhbmdsZXMuaW5zZXJ0KHBvaW50c1tpXSwgdXNlUmFuZG9tKVxuICB9XG4gIFxuICAvL0V4dHJhY3QgYm91bmRhcnkgY2VsbHNcbiAgcmV0dXJuIHRyaWFuZ2xlcy5ib3VuZGFyeSgpXG59IiwiXCJ1c2Ugc3RyaWN0XCJcblxubW9kdWxlLmV4cG9ydHMgPSBmYXN0VHdvU3VtXG5cbmZ1bmN0aW9uIGZhc3RUd29TdW0oYSwgYiwgcmVzdWx0KSB7XG5cdHZhciB4ID0gYSArIGJcblx0dmFyIGJ2ID0geCAtIGFcblx0dmFyIGF2ID0geCAtIGJ2XG5cdHZhciBiciA9IGIgLSBidlxuXHR2YXIgYXIgPSBhIC0gYXZcblx0aWYocmVzdWx0KSB7XG5cdFx0cmVzdWx0WzBdID0gYXIgKyBiclxuXHRcdHJlc3VsdFsxXSA9IHhcblx0XHRyZXR1cm4gcmVzdWx0XG5cdH1cblx0cmV0dXJuIFthciticiwgeF1cbn0iLCJcInVzZSBzdHJpY3RcIlxuXG52YXIgdHdvUHJvZHVjdCA9IHJlcXVpcmUoXCJ0d28tcHJvZHVjdFwiKVxudmFyIHR3b1N1bSA9IHJlcXVpcmUoXCJ0d28tc3VtXCIpXG5cbm1vZHVsZS5leHBvcnRzID0gc2NhbGVMaW5lYXJFeHBhbnNpb25cblxuZnVuY3Rpb24gc2NhbGVMaW5lYXJFeHBhbnNpb24oZSwgc2NhbGUpIHtcbiAgdmFyIG4gPSBlLmxlbmd0aFxuICBpZihuID09PSAxKSB7XG4gICAgdmFyIHRzID0gdHdvUHJvZHVjdChlWzBdLCBzY2FsZSlcbiAgICBpZih0c1swXSkge1xuICAgICAgcmV0dXJuIHRzXG4gICAgfVxuICAgIHJldHVybiBbIHRzWzFdIF1cbiAgfVxuICB2YXIgZyA9IG5ldyBBcnJheSgyICogbilcbiAgdmFyIHEgPSBbMC4xLCAwLjFdXG4gIHZhciB0ID0gWzAuMSwgMC4xXVxuICB2YXIgY291bnQgPSAwXG4gIHR3b1Byb2R1Y3QoZVswXSwgc2NhbGUsIHEpXG4gIGlmKHFbMF0pIHtcbiAgICBnW2NvdW50KytdID0gcVswXVxuICB9XG4gIGZvcih2YXIgaT0xOyBpPG47ICsraSkge1xuICAgIHR3b1Byb2R1Y3QoZVtpXSwgc2NhbGUsIHQpXG4gICAgdmFyIHBxID0gcVsxXVxuICAgIHR3b1N1bShwcSwgdFswXSwgcSlcbiAgICBpZihxWzBdKSB7XG4gICAgICBnW2NvdW50KytdID0gcVswXVxuICAgIH1cbiAgICB2YXIgYSA9IHRbMV1cbiAgICB2YXIgYiA9IHFbMV1cbiAgICB2YXIgeCA9IGEgKyBiXG4gICAgdmFyIGJ2ID0geCAtIGFcbiAgICB2YXIgeSA9IGIgLSBidlxuICAgIHFbMV0gPSB4XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gIH1cbiAgaWYocVsxXSkge1xuICAgIGdbY291bnQrK10gPSBxWzFdXG4gIH1cbiAgaWYoY291bnQgPT09IDApIHtcbiAgICBnW2NvdW50KytdID0gMC4wXG4gIH1cbiAgZy5sZW5ndGggPSBjb3VudFxuICByZXR1cm4gZ1xufSIsIlwidXNlIHN0cmljdFwiXG5cbm1vZHVsZS5leHBvcnRzID0gcm9idXN0U3VidHJhY3RcblxuLy9FYXN5IGNhc2U6IEFkZCB0d28gc2NhbGFyc1xuZnVuY3Rpb24gc2NhbGFyU2NhbGFyKGEsIGIpIHtcbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgYXYgPSB4IC0gYnZcbiAgdmFyIGJyID0gYiAtIGJ2XG4gIHZhciBhciA9IGEgLSBhdlxuICB2YXIgeSA9IGFyICsgYnJcbiAgaWYoeSkge1xuICAgIHJldHVybiBbeSwgeF1cbiAgfVxuICByZXR1cm4gW3hdXG59XG5cbmZ1bmN0aW9uIHJvYnVzdFN1YnRyYWN0KGUsIGYpIHtcbiAgdmFyIG5lID0gZS5sZW5ndGh8MFxuICB2YXIgbmYgPSBmLmxlbmd0aHwwXG4gIGlmKG5lID09PSAxICYmIG5mID09PSAxKSB7XG4gICAgcmV0dXJuIHNjYWxhclNjYWxhcihlWzBdLCAtZlswXSlcbiAgfVxuICB2YXIgbiA9IG5lICsgbmZcbiAgdmFyIGcgPSBuZXcgQXJyYXkobilcbiAgdmFyIGNvdW50ID0gMFxuICB2YXIgZXB0ciA9IDBcbiAgdmFyIGZwdHIgPSAwXG4gIHZhciBhYnMgPSBNYXRoLmFic1xuICB2YXIgZWkgPSBlW2VwdHJdXG4gIHZhciBlYSA9IGFicyhlaSlcbiAgdmFyIGZpID0gLWZbZnB0cl1cbiAgdmFyIGZhID0gYWJzKGZpKVxuICB2YXIgYSwgYlxuICBpZihlYSA8IGZhKSB7XG4gICAgYiA9IGVpXG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICAgIGVhID0gYWJzKGVpKVxuICAgIH1cbiAgfSBlbHNlIHtcbiAgICBiID0gZmlcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gLWZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICBpZigoZXB0ciA8IG5lICYmIGVhIDwgZmEpIHx8IChmcHRyID49IG5mKSkge1xuICAgIGEgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYSA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IC1mW2ZwdHJdXG4gICAgICBmYSA9IGFicyhmaSlcbiAgICB9XG4gIH1cbiAgdmFyIHggPSBhICsgYlxuICB2YXIgYnYgPSB4IC0gYVxuICB2YXIgeSA9IGIgLSBidlxuICB2YXIgcTAgPSB5XG4gIHZhciBxMSA9IHhcbiAgdmFyIF94LCBfYnYsIF9hdiwgX2JyLCBfYXJcbiAgd2hpbGUoZXB0ciA8IG5lICYmIGZwdHIgPCBuZikge1xuICAgIGlmKGVhIDwgZmEpIHtcbiAgICAgIGEgPSBlaVxuICAgICAgZXB0ciArPSAxXG4gICAgICBpZihlcHRyIDwgbmUpIHtcbiAgICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICAgIGVhID0gYWJzKGVpKVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICBhID0gZmlcbiAgICAgIGZwdHIgKz0gMVxuICAgICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICAgIGZpID0gLWZbZnB0cl1cbiAgICAgICAgZmEgPSBhYnMoZmkpXG4gICAgICB9XG4gICAgfVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9XG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgfVxuICB3aGlsZShlcHRyIDwgbmUpIHtcbiAgICBhID0gZWlcbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZXB0ciArPSAxXG4gICAgaWYoZXB0ciA8IG5lKSB7XG4gICAgICBlaSA9IGVbZXB0cl1cbiAgICB9XG4gIH1cbiAgd2hpbGUoZnB0ciA8IG5mKSB7XG4gICAgYSA9IGZpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH0gXG4gICAgX3ggPSBxMSArIHhcbiAgICBfYnYgPSBfeCAtIHExXG4gICAgX2F2ID0gX3ggLSBfYnZcbiAgICBfYnIgPSB4IC0gX2J2XG4gICAgX2FyID0gcTEgLSBfYXZcbiAgICBxMCA9IF9hciArIF9iclxuICAgIHExID0gX3hcbiAgICBmcHRyICs9IDFcbiAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgIGZpID0gLWZbZnB0cl1cbiAgICB9XG4gIH1cbiAgaWYocTApIHtcbiAgICBnW2NvdW50KytdID0gcTBcbiAgfVxuICBpZihxMSkge1xuICAgIGdbY291bnQrK10gPSBxMVxuICB9XG4gIGlmKCFjb3VudCkge1xuICAgIGdbY291bnQrK10gPSAwLjAgIFxuICB9XG4gIGcubGVuZ3RoID0gY291bnRcbiAgcmV0dXJuIGdcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IGxpbmVhckV4cGFuc2lvblN1bVxuXG4vL0Vhc3kgY2FzZTogQWRkIHR3byBzY2FsYXJzXG5mdW5jdGlvbiBzY2FsYXJTY2FsYXIoYSwgYikge1xuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciBhdiA9IHggLSBidlxuICB2YXIgYnIgPSBiIC0gYnZcbiAgdmFyIGFyID0gYSAtIGF2XG4gIHZhciB5ID0gYXIgKyBiclxuICBpZih5KSB7XG4gICAgcmV0dXJuIFt5LCB4XVxuICB9XG4gIHJldHVybiBbeF1cbn1cblxuZnVuY3Rpb24gbGluZWFyRXhwYW5zaW9uU3VtKGUsIGYpIHtcbiAgdmFyIG5lID0gZS5sZW5ndGh8MFxuICB2YXIgbmYgPSBmLmxlbmd0aHwwXG4gIGlmKG5lID09PSAxICYmIG5mID09PSAxKSB7XG4gICAgcmV0dXJuIHNjYWxhclNjYWxhcihlWzBdLCBmWzBdKVxuICB9XG4gIHZhciBuID0gbmUgKyBuZlxuICB2YXIgZyA9IG5ldyBBcnJheShuKVxuICB2YXIgY291bnQgPSAwXG4gIHZhciBlcHRyID0gMFxuICB2YXIgZnB0ciA9IDBcbiAgdmFyIGFicyA9IE1hdGguYWJzXG4gIHZhciBlaSA9IGVbZXB0cl1cbiAgdmFyIGVhID0gYWJzKGVpKVxuICB2YXIgZmkgPSBmW2ZwdHJdXG4gIHZhciBmYSA9IGFicyhmaSlcbiAgdmFyIGEsIGJcbiAgaWYoZWEgPCBmYSkge1xuICAgIGIgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYiA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IGZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICBpZigoZXB0ciA8IG5lICYmIGVhIDwgZmEpIHx8IChmcHRyID49IG5mKSkge1xuICAgIGEgPSBlaVxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgICBlYSA9IGFicyhlaSlcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgYSA9IGZpXG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IGZbZnB0cl1cbiAgICAgIGZhID0gYWJzKGZpKVxuICAgIH1cbiAgfVxuICB2YXIgeCA9IGEgKyBiXG4gIHZhciBidiA9IHggLSBhXG4gIHZhciB5ID0gYiAtIGJ2XG4gIHZhciBxMCA9IHlcbiAgdmFyIHExID0geFxuICB2YXIgX3gsIF9idiwgX2F2LCBfYnIsIF9hclxuICB3aGlsZShlcHRyIDwgbmUgJiYgZnB0ciA8IG5mKSB7XG4gICAgaWYoZWEgPCBmYSkge1xuICAgICAgYSA9IGVpXG4gICAgICBlcHRyICs9IDFcbiAgICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgICBlaSA9IGVbZXB0cl1cbiAgICAgICAgZWEgPSBhYnMoZWkpXG4gICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgIGEgPSBmaVxuICAgICAgZnB0ciArPSAxXG4gICAgICBpZihmcHRyIDwgbmYpIHtcbiAgICAgICAgZmkgPSBmW2ZwdHJdXG4gICAgICAgIGZhID0gYWJzKGZpKVxuICAgICAgfVxuICAgIH1cbiAgICBiID0gcTBcbiAgICB4ID0gYSArIGJcbiAgICBidiA9IHggLSBhXG4gICAgeSA9IGIgLSBidlxuICAgIGlmKHkpIHtcbiAgICAgIGdbY291bnQrK10gPSB5XG4gICAgfVxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gIH1cbiAgd2hpbGUoZXB0ciA8IG5lKSB7XG4gICAgYSA9IGVpXG4gICAgYiA9IHEwXG4gICAgeCA9IGEgKyBiXG4gICAgYnYgPSB4IC0gYVxuICAgIHkgPSBiIC0gYnZcbiAgICBpZih5KSB7XG4gICAgICBnW2NvdW50KytdID0geVxuICAgIH1cbiAgICBfeCA9IHExICsgeFxuICAgIF9idiA9IF94IC0gcTFcbiAgICBfYXYgPSBfeCAtIF9idlxuICAgIF9iciA9IHggLSBfYnZcbiAgICBfYXIgPSBxMSAtIF9hdlxuICAgIHEwID0gX2FyICsgX2JyXG4gICAgcTEgPSBfeFxuICAgIGVwdHIgKz0gMVxuICAgIGlmKGVwdHIgPCBuZSkge1xuICAgICAgZWkgPSBlW2VwdHJdXG4gICAgfVxuICB9XG4gIHdoaWxlKGZwdHIgPCBuZikge1xuICAgIGEgPSBmaVxuICAgIGIgPSBxMFxuICAgIHggPSBhICsgYlxuICAgIGJ2ID0geCAtIGFcbiAgICB5ID0gYiAtIGJ2XG4gICAgaWYoeSkge1xuICAgICAgZ1tjb3VudCsrXSA9IHlcbiAgICB9IFxuICAgIF94ID0gcTEgKyB4XG4gICAgX2J2ID0gX3ggLSBxMVxuICAgIF9hdiA9IF94IC0gX2J2XG4gICAgX2JyID0geCAtIF9idlxuICAgIF9hciA9IHExIC0gX2F2XG4gICAgcTAgPSBfYXIgKyBfYnJcbiAgICBxMSA9IF94XG4gICAgZnB0ciArPSAxXG4gICAgaWYoZnB0ciA8IG5mKSB7XG4gICAgICBmaSA9IGZbZnB0cl1cbiAgICB9XG4gIH1cbiAgaWYocTApIHtcbiAgICBnW2NvdW50KytdID0gcTBcbiAgfVxuICBpZihxMSkge1xuICAgIGdbY291bnQrK10gPSBxMVxuICB9XG4gIGlmKCFjb3VudCkge1xuICAgIGdbY291bnQrK10gPSAwLjAgIFxuICB9XG4gIGcubGVuZ3RoID0gY291bnRcbiAgcmV0dXJuIGdcbn0iLCJcInVzZSBzdHJpY3RcIlxuXG5tb2R1bGUuZXhwb3J0cyA9IHR3b1Byb2R1Y3RcblxudmFyIFNQTElUVEVSID0gKyhNYXRoLnBvdygyLCAyNykgKyAxLjApXG5cbmZ1bmN0aW9uIHR3b1Byb2R1Y3QoYSwgYiwgcmVzdWx0KSB7XG4gIHZhciB4ID0gYSAqIGJcblxuICB2YXIgYyA9IFNQTElUVEVSICogYVxuICB2YXIgYWJpZyA9IGMgLSBhXG4gIHZhciBhaGkgPSBjIC0gYWJpZ1xuICB2YXIgYWxvID0gYSAtIGFoaVxuXG4gIHZhciBkID0gU1BMSVRURVIgKiBiXG4gIHZhciBiYmlnID0gZCAtIGJcbiAgdmFyIGJoaSA9IGQgLSBiYmlnXG4gIHZhciBibG8gPSBiIC0gYmhpXG5cbiAgdmFyIGVycjEgPSB4IC0gKGFoaSAqIGJoaSlcbiAgdmFyIGVycjIgPSBlcnIxIC0gKGFsbyAqIGJoaSlcbiAgdmFyIGVycjMgPSBlcnIyIC0gKGFoaSAqIGJsbylcblxuICB2YXIgeSA9IGFsbyAqIGJsbyAtIGVycjNcblxuICBpZihyZXN1bHQpIHtcbiAgICByZXN1bHRbMF0gPSB5XG4gICAgcmVzdWx0WzFdID0geFxuICAgIHJldHVybiByZXN1bHRcbiAgfVxuXG4gIHJldHVybiBbIHksIHggXVxufSIsIlwidXNlIHN0cmljdFwiXG5cbnZhciB0d29Qcm9kdWN0ID0gcmVxdWlyZShcInR3by1wcm9kdWN0XCIpXG52YXIgcm9idXN0U3VtID0gcmVxdWlyZShcInJvYnVzdC1zdW1cIilcbnZhciByb2J1c3RTY2FsZSA9IHJlcXVpcmUoXCJyb2J1c3Qtc2NhbGVcIilcbnZhciByb2J1c3RTdWJ0cmFjdCA9IHJlcXVpcmUoXCJyb2J1c3Qtc3VidHJhY3RcIilcblxudmFyIE5VTV9FWFBBTkQgPSA1XG5cbnZhciBFUFNJTE9OICAgICA9IDEuMTEwMjIzMDI0NjI1MTU2NWUtMTZcbnZhciBFUlJCT1VORDMgICA9ICgzLjAgKyAxNi4wICogRVBTSUxPTikgKiBFUFNJTE9OXG52YXIgRVJSQk9VTkQ0ICAgPSAoNy4wICsgNTYuMCAqIEVQU0lMT04pICogRVBTSUxPTlxuXG5mdW5jdGlvbiBjb2ZhY3RvcihtLCBjKSB7XG4gIHZhciByZXN1bHQgPSBuZXcgQXJyYXkobS5sZW5ndGgtMSlcbiAgZm9yKHZhciBpPTE7IGk8bS5sZW5ndGg7ICsraSkge1xuICAgIHZhciByID0gcmVzdWx0W2ktMV0gPSBuZXcgQXJyYXkobS5sZW5ndGgtMSlcbiAgICBmb3IodmFyIGo9MCxrPTA7IGo8bS5sZW5ndGg7ICsraikge1xuICAgICAgaWYoaiA9PT0gYykge1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgcltrKytdID0gbVtpXVtqXVxuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzdWx0XG59XG5cbmZ1bmN0aW9uIG1hdHJpeChuKSB7XG4gIHZhciByZXN1bHQgPSBuZXcgQXJyYXkobilcbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgcmVzdWx0W2ldID0gbmV3IEFycmF5KG4pXG4gICAgZm9yKHZhciBqPTA7IGo8bjsgKytqKSB7XG4gICAgICByZXN1bHRbaV1bal0gPSBbXCJtXCIsIGosIFwiW1wiLCAobi1pLTEpLCBcIl1cIl0uam9pbihcIlwiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzdWx0XG59XG5cbmZ1bmN0aW9uIHNpZ24obikge1xuICBpZihuICYgMSkge1xuICAgIHJldHVybiBcIi1cIlxuICB9XG4gIHJldHVybiBcIlwiXG59XG5cbmZ1bmN0aW9uIGdlbmVyYXRlU3VtKGV4cHIpIHtcbiAgaWYoZXhwci5sZW5ndGggPT09IDEpIHtcbiAgICByZXR1cm4gZXhwclswXVxuICB9IGVsc2UgaWYoZXhwci5sZW5ndGggPT09IDIpIHtcbiAgICByZXR1cm4gW1wic3VtKFwiLCBleHByWzBdLCBcIixcIiwgZXhwclsxXSwgXCIpXCJdLmpvaW4oXCJcIilcbiAgfSBlbHNlIHtcbiAgICB2YXIgbSA9IGV4cHIubGVuZ3RoPj4xXG4gICAgcmV0dXJuIFtcInN1bShcIiwgZ2VuZXJhdGVTdW0oZXhwci5zbGljZSgwLCBtKSksIFwiLFwiLCBnZW5lcmF0ZVN1bShleHByLnNsaWNlKG0pKSwgXCIpXCJdLmpvaW4oXCJcIilcbiAgfVxufVxuXG5mdW5jdGlvbiBkZXRlcm1pbmFudChtKSB7XG4gIGlmKG0ubGVuZ3RoID09PSAyKSB7XG4gICAgcmV0dXJuIFtbXCJzdW0ocHJvZChcIiwgbVswXVswXSwgXCIsXCIsIG1bMV1bMV0sIFwiKSxwcm9kKC1cIiwgbVswXVsxXSwgXCIsXCIsIG1bMV1bMF0sIFwiKSlcIl0uam9pbihcIlwiKV1cbiAgfSBlbHNlIHtcbiAgICB2YXIgZXhwciA9IFtdXG4gICAgZm9yKHZhciBpPTA7IGk8bS5sZW5ndGg7ICsraSkge1xuICAgICAgZXhwci5wdXNoKFtcInNjYWxlKFwiLCBnZW5lcmF0ZVN1bShkZXRlcm1pbmFudChjb2ZhY3RvcihtLCBpKSkpLCBcIixcIiwgc2lnbihpKSwgbVswXVtpXSwgXCIpXCJdLmpvaW4oXCJcIikpXG4gICAgfVxuICAgIHJldHVybiBleHByXG4gIH1cbn1cblxuZnVuY3Rpb24gb3JpZW50YXRpb24obikge1xuICB2YXIgcG9zID0gW11cbiAgdmFyIG5lZyA9IFtdXG4gIHZhciBtID0gbWF0cml4KG4pXG4gIHZhciBhcmdzID0gW11cbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgaWYoKGkmMSk9PT0wKSB7XG4gICAgICBwb3MucHVzaC5hcHBseShwb3MsIGRldGVybWluYW50KGNvZmFjdG9yKG0sIGkpKSlcbiAgICB9IGVsc2Uge1xuICAgICAgbmVnLnB1c2guYXBwbHkobmVnLCBkZXRlcm1pbmFudChjb2ZhY3RvcihtLCBpKSkpXG4gICAgfVxuICAgIGFyZ3MucHVzaChcIm1cIiArIGkpXG4gIH1cbiAgdmFyIHBvc0V4cHIgPSBnZW5lcmF0ZVN1bShwb3MpXG4gIHZhciBuZWdFeHByID0gZ2VuZXJhdGVTdW0obmVnKVxuICB2YXIgZnVuY05hbWUgPSBcIm9yaWVudGF0aW9uXCIgKyBuICsgXCJFeGFjdFwiXG4gIHZhciBjb2RlID0gW1wiZnVuY3Rpb24gXCIsIGZ1bmNOYW1lLCBcIihcIiwgYXJncy5qb2luKCksIFwiKXt2YXIgcD1cIiwgcG9zRXhwciwgXCIsbj1cIiwgbmVnRXhwciwgXCIsZD1zdWIocCxuKTtcXFxucmV0dXJuIGRbZC5sZW5ndGgtMV07fTtyZXR1cm4gXCIsIGZ1bmNOYW1lXS5qb2luKFwiXCIpXG4gIHZhciBwcm9jID0gbmV3IEZ1bmN0aW9uKFwic3VtXCIsIFwicHJvZFwiLCBcInNjYWxlXCIsIFwic3ViXCIsIGNvZGUpXG4gIHJldHVybiBwcm9jKHJvYnVzdFN1bSwgdHdvUHJvZHVjdCwgcm9idXN0U2NhbGUsIHJvYnVzdFN1YnRyYWN0KVxufVxuXG52YXIgb3JpZW50YXRpb24zRXhhY3QgPSBvcmllbnRhdGlvbigzKVxudmFyIG9yaWVudGF0aW9uNEV4YWN0ID0gb3JpZW50YXRpb24oNClcblxudmFyIENBQ0hFRCA9IFtcbiAgZnVuY3Rpb24gb3JpZW50YXRpb24wKCkgeyByZXR1cm4gMCB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjEoKSB7IHJldHVybiAwIH0sXG4gIGZ1bmN0aW9uIG9yaWVudGF0aW9uMihhLCBiKSB7IFxuICAgIHJldHVybiBiWzBdIC0gYVswXVxuICB9LFxuICBmdW5jdGlvbiBvcmllbnRhdGlvbjMoYSwgYiwgYykge1xuICAgIHZhciBsID0gKGFbMV0gLSBjWzFdKSAqIChiWzBdIC0gY1swXSlcbiAgICB2YXIgciA9IChhWzBdIC0gY1swXSkgKiAoYlsxXSAtIGNbMV0pXG4gICAgdmFyIGRldCA9IGwgLSByXG4gICAgdmFyIHNcbiAgICBpZihsID4gMCkge1xuICAgICAgaWYociA8PSAwKSB7XG4gICAgICAgIHJldHVybiBkZXRcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIHMgPSBsICsgclxuICAgICAgfVxuICAgIH0gZWxzZSBpZihsIDwgMCkge1xuICAgICAgaWYociA+PSAwKSB7XG4gICAgICAgIHJldHVybiBkZXRcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIHMgPSAtKGwgKyByKVxuICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICByZXR1cm4gZGV0XG4gICAgfVxuICAgIHZhciB0b2wgPSBFUlJCT1VORDMgKiBzXG4gICAgaWYoZGV0ID49IHRvbCB8fCBkZXQgPD0gLXRvbCkge1xuICAgICAgcmV0dXJuIGRldFxuICAgIH1cbiAgICByZXR1cm4gb3JpZW50YXRpb24zRXhhY3QoYSwgYiwgYylcbiAgfSxcbiAgZnVuY3Rpb24gb3JpZW50YXRpb240KGEsYixjLGQpIHtcbiAgICB2YXIgYWR4ID0gYVswXSAtIGRbMF1cbiAgICB2YXIgYmR4ID0gYlswXSAtIGRbMF1cbiAgICB2YXIgY2R4ID0gY1swXSAtIGRbMF1cbiAgICB2YXIgYWR5ID0gYVsxXSAtIGRbMV1cbiAgICB2YXIgYmR5ID0gYlsxXSAtIGRbMV1cbiAgICB2YXIgY2R5ID0gY1sxXSAtIGRbMV1cbiAgICB2YXIgYWR6ID0gYVsyXSAtIGRbMl1cbiAgICB2YXIgYmR6ID0gYlsyXSAtIGRbMl1cbiAgICB2YXIgY2R6ID0gY1syXSAtIGRbMl1cbiAgICB2YXIgYmR4Y2R5ID0gYmR4ICogY2R5XG4gICAgdmFyIGNkeGJkeSA9IGNkeCAqIGJkeVxuICAgIHZhciBjZHhhZHkgPSBjZHggKiBhZHlcbiAgICB2YXIgYWR4Y2R5ID0gYWR4ICogY2R5XG4gICAgdmFyIGFkeGJkeSA9IGFkeCAqIGJkeVxuICAgIHZhciBiZHhhZHkgPSBiZHggKiBhZHlcbiAgICB2YXIgZGV0ID0gYWR6ICogKGJkeGNkeSAtIGNkeGJkeSkgXG4gICAgICAgICAgICArIGJkeiAqIChjZHhhZHkgLSBhZHhjZHkpXG4gICAgICAgICAgICArIGNkeiAqIChhZHhiZHkgLSBiZHhhZHkpXG4gICAgdmFyIHBlcm1hbmVudCA9IChNYXRoLmFicyhiZHhjZHkpICsgTWF0aC5hYnMoY2R4YmR5KSkgKiBNYXRoLmFicyhhZHopXG4gICAgICAgICAgICAgICAgICArIChNYXRoLmFicyhjZHhhZHkpICsgTWF0aC5hYnMoYWR4Y2R5KSkgKiBNYXRoLmFicyhiZHopXG4gICAgICAgICAgICAgICAgICArIChNYXRoLmFicyhhZHhiZHkpICsgTWF0aC5hYnMoYmR4YWR5KSkgKiBNYXRoLmFicyhjZHopXG4gICAgdmFyIHRvbCA9IEVSUkJPVU5ENCAqIHBlcm1hbmVudFxuICAgIGlmICgoZGV0ID4gdG9sKSB8fCAoLWRldCA+IHRvbCkpIHtcbiAgICAgIHJldHVybiBkZXRcbiAgICB9XG4gICAgcmV0dXJuIG9yaWVudGF0aW9uNEV4YWN0KGEsYixjLGQpXG4gIH1cbl1cblxuZnVuY3Rpb24gc2xvd09yaWVudChhcmdzKSB7XG4gIHZhciBwcm9jID0gQ0FDSEVEW2FyZ3MubGVuZ3RoXVxuICBpZighcHJvYykge1xuICAgIHByb2MgPSBDQUNIRURbYXJncy5sZW5ndGhdID0gb3JpZW50YXRpb24oYXJncy5sZW5ndGgpXG4gIH1cbiAgcmV0dXJuIHByb2MuYXBwbHkodW5kZWZpbmVkLCBhcmdzKVxufVxuXG5mdW5jdGlvbiBnZW5lcmF0ZU9yaWVudGF0aW9uUHJvYygpIHtcbiAgd2hpbGUoQ0FDSEVELmxlbmd0aCA8PSBOVU1fRVhQQU5EKSB7XG4gICAgQ0FDSEVELnB1c2gob3JpZW50YXRpb24oQ0FDSEVELmxlbmd0aCkpXG4gIH1cbiAgdmFyIGFyZ3MgPSBbXVxuICB2YXIgcHJvY0FyZ3MgPSBbXCJzbG93XCJdXG4gIGZvcih2YXIgaT0wOyBpPD1OVU1fRVhQQU5EOyArK2kpIHtcbiAgICBhcmdzLnB1c2goXCJhXCIgKyBpKVxuICAgIHByb2NBcmdzLnB1c2goXCJvXCIgKyBpKVxuICB9XG4gIHZhciBjb2RlID0gW1xuICAgIFwiZnVuY3Rpb24gZ2V0T3JpZW50YXRpb24oXCIsIGFyZ3Muam9pbigpLCBcIil7c3dpdGNoKGFyZ3VtZW50cy5sZW5ndGgpe2Nhc2UgMDpjYXNlIDE6cmV0dXJuIDA7XCJcbiAgXVxuICBmb3IodmFyIGk9MjsgaTw9TlVNX0VYUEFORDsgKytpKSB7XG4gICAgY29kZS5wdXNoKFwiY2FzZSBcIiwgaSwgXCI6cmV0dXJuIG9cIiwgaSwgXCIoXCIsIGFyZ3Muc2xpY2UoMCwgaSkuam9pbigpLCBcIik7XCIpXG4gIH1cbiAgY29kZS5wdXNoKFwifXZhciBzPW5ldyBBcnJheShhcmd1bWVudHMubGVuZ3RoKTtmb3IodmFyIGk9MDtpPGFyZ3VtZW50cy5sZW5ndGg7KytpKXtzW2ldPWFyZ3VtZW50c1tpXX07cmV0dXJuIHNsb3cocyk7fXJldHVybiBnZXRPcmllbnRhdGlvblwiKVxuICBwcm9jQXJncy5wdXNoKGNvZGUuam9pbihcIlwiKSlcblxuICB2YXIgcHJvYyA9IEZ1bmN0aW9uLmFwcGx5KHVuZGVmaW5lZCwgcHJvY0FyZ3MpXG4gIG1vZHVsZS5leHBvcnRzID0gcHJvYy5hcHBseSh1bmRlZmluZWQsIFtzbG93T3JpZW50XS5jb25jYXQoQ0FDSEVEKSlcbiAgZm9yKHZhciBpPTA7IGk8PU5VTV9FWFBBTkQ7ICsraSkge1xuICAgIG1vZHVsZS5leHBvcnRzW2ldID0gQ0FDSEVEW2ldXG4gIH1cbn1cblxuZ2VuZXJhdGVPcmllbnRhdGlvblByb2MoKSIsIi8qKlxuICogQml0IHR3aWRkbGluZyBoYWNrcyBmb3IgSmF2YVNjcmlwdC5cbiAqXG4gKiBBdXRob3I6IE1pa29sYSBMeXNlbmtvXG4gKlxuICogUG9ydGVkIGZyb20gU3RhbmZvcmQgYml0IHR3aWRkbGluZyBoYWNrIGxpYnJhcnk6XG4gKiAgICBodHRwOi8vZ3JhcGhpY3Muc3RhbmZvcmQuZWR1L35zZWFuZGVyL2JpdGhhY2tzLmh0bWxcbiAqL1xuXG5cInVzZSBzdHJpY3RcIjsgXCJ1c2UgcmVzdHJpY3RcIjtcblxuLy9OdW1iZXIgb2YgYml0cyBpbiBhbiBpbnRlZ2VyXG52YXIgSU5UX0JJVFMgPSAzMjtcblxuLy9Db25zdGFudHNcbmV4cG9ydHMuSU5UX0JJVFMgID0gSU5UX0JJVFM7XG5leHBvcnRzLklOVF9NQVggICA9ICAweDdmZmZmZmZmO1xuZXhwb3J0cy5JTlRfTUlOICAgPSAtMTw8KElOVF9CSVRTLTEpO1xuXG4vL1JldHVybnMgLTEsIDAsICsxIGRlcGVuZGluZyBvbiBzaWduIG9mIHhcbmV4cG9ydHMuc2lnbiA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICh2ID4gMCkgLSAodiA8IDApO1xufVxuXG4vL0NvbXB1dGVzIGFic29sdXRlIHZhbHVlIG9mIGludGVnZXJcbmV4cG9ydHMuYWJzID0gZnVuY3Rpb24odikge1xuICB2YXIgbWFzayA9IHYgPj4gKElOVF9CSVRTLTEpO1xuICByZXR1cm4gKHYgXiBtYXNrKSAtIG1hc2s7XG59XG5cbi8vQ29tcHV0ZXMgbWluaW11bSBvZiBpbnRlZ2VycyB4IGFuZCB5XG5leHBvcnRzLm1pbiA9IGZ1bmN0aW9uKHgsIHkpIHtcbiAgcmV0dXJuIHkgXiAoKHggXiB5KSAmIC0oeCA8IHkpKTtcbn1cblxuLy9Db21wdXRlcyBtYXhpbXVtIG9mIGludGVnZXJzIHggYW5kIHlcbmV4cG9ydHMubWF4ID0gZnVuY3Rpb24oeCwgeSkge1xuICByZXR1cm4geCBeICgoeCBeIHkpICYgLSh4IDwgeSkpO1xufVxuXG4vL0NoZWNrcyBpZiBhIG51bWJlciBpcyBhIHBvd2VyIG9mIHR3b1xuZXhwb3J0cy5pc1BvdzIgPSBmdW5jdGlvbih2KSB7XG4gIHJldHVybiAhKHYgJiAodi0xKSkgJiYgKCEhdik7XG59XG5cbi8vQ29tcHV0ZXMgbG9nIGJhc2UgMiBvZiB2XG5leHBvcnRzLmxvZzIgPSBmdW5jdGlvbih2KSB7XG4gIHZhciByLCBzaGlmdDtcbiAgciA9ICAgICAodiA+IDB4RkZGRikgPDwgNDsgdiA+Pj49IHI7XG4gIHNoaWZ0ID0gKHYgPiAweEZGICApIDw8IDM7IHYgPj4+PSBzaGlmdDsgciB8PSBzaGlmdDtcbiAgc2hpZnQgPSAodiA+IDB4RiAgICkgPDwgMjsgdiA+Pj49IHNoaWZ0OyByIHw9IHNoaWZ0O1xuICBzaGlmdCA9ICh2ID4gMHgzICAgKSA8PCAxOyB2ID4+Pj0gc2hpZnQ7IHIgfD0gc2hpZnQ7XG4gIHJldHVybiByIHwgKHYgPj4gMSk7XG59XG5cbi8vQ29tcHV0ZXMgbG9nIGJhc2UgMTAgb2YgdlxuZXhwb3J0cy5sb2cxMCA9IGZ1bmN0aW9uKHYpIHtcbiAgcmV0dXJuICAodiA+PSAxMDAwMDAwMDAwKSA/IDkgOiAodiA+PSAxMDAwMDAwMDApID8gOCA6ICh2ID49IDEwMDAwMDAwKSA/IDcgOlxuICAgICAgICAgICh2ID49IDEwMDAwMDApID8gNiA6ICh2ID49IDEwMDAwMCkgPyA1IDogKHYgPj0gMTAwMDApID8gNCA6XG4gICAgICAgICAgKHYgPj0gMTAwMCkgPyAzIDogKHYgPj0gMTAwKSA/IDIgOiAodiA+PSAxMCkgPyAxIDogMDtcbn1cblxuLy9Db3VudHMgbnVtYmVyIG9mIGJpdHNcbmV4cG9ydHMucG9wQ291bnQgPSBmdW5jdGlvbih2KSB7XG4gIHYgPSB2IC0gKCh2ID4+PiAxKSAmIDB4NTU1NTU1NTUpO1xuICB2ID0gKHYgJiAweDMzMzMzMzMzKSArICgodiA+Pj4gMikgJiAweDMzMzMzMzMzKTtcbiAgcmV0dXJuICgodiArICh2ID4+PiA0KSAmIDB4RjBGMEYwRikgKiAweDEwMTAxMDEpID4+PiAyNDtcbn1cblxuLy9Db3VudHMgbnVtYmVyIG9mIHRyYWlsaW5nIHplcm9zXG5mdW5jdGlvbiBjb3VudFRyYWlsaW5nWmVyb3Modikge1xuICB2YXIgYyA9IDMyO1xuICB2ICY9IC12O1xuICBpZiAodikgYy0tO1xuICBpZiAodiAmIDB4MDAwMEZGRkYpIGMgLT0gMTY7XG4gIGlmICh2ICYgMHgwMEZGMDBGRikgYyAtPSA4O1xuICBpZiAodiAmIDB4MEYwRjBGMEYpIGMgLT0gNDtcbiAgaWYgKHYgJiAweDMzMzMzMzMzKSBjIC09IDI7XG4gIGlmICh2ICYgMHg1NTU1NTU1NSkgYyAtPSAxO1xuICByZXR1cm4gYztcbn1cbmV4cG9ydHMuY291bnRUcmFpbGluZ1plcm9zID0gY291bnRUcmFpbGluZ1plcm9zO1xuXG4vL1JvdW5kcyB0byBuZXh0IHBvd2VyIG9mIDJcbmV4cG9ydHMubmV4dFBvdzIgPSBmdW5jdGlvbih2KSB7XG4gIHYgKz0gdiA9PT0gMDtcbiAgLS12O1xuICB2IHw9IHYgPj4+IDE7XG4gIHYgfD0gdiA+Pj4gMjtcbiAgdiB8PSB2ID4+PiA0O1xuICB2IHw9IHYgPj4+IDg7XG4gIHYgfD0gdiA+Pj4gMTY7XG4gIHJldHVybiB2ICsgMTtcbn1cblxuLy9Sb3VuZHMgZG93biB0byBwcmV2aW91cyBwb3dlciBvZiAyXG5leHBvcnRzLnByZXZQb3cyID0gZnVuY3Rpb24odikge1xuICB2IHw9IHYgPj4+IDE7XG4gIHYgfD0gdiA+Pj4gMjtcbiAgdiB8PSB2ID4+PiA0O1xuICB2IHw9IHYgPj4+IDg7XG4gIHYgfD0gdiA+Pj4gMTY7XG4gIHJldHVybiB2IC0gKHY+Pj4xKTtcbn1cblxuLy9Db21wdXRlcyBwYXJpdHkgb2Ygd29yZFxuZXhwb3J0cy5wYXJpdHkgPSBmdW5jdGlvbih2KSB7XG4gIHYgXj0gdiA+Pj4gMTY7XG4gIHYgXj0gdiA+Pj4gODtcbiAgdiBePSB2ID4+PiA0O1xuICB2ICY9IDB4ZjtcbiAgcmV0dXJuICgweDY5OTYgPj4+IHYpICYgMTtcbn1cblxudmFyIFJFVkVSU0VfVEFCTEUgPSBuZXcgQXJyYXkoMjU2KTtcblxuKGZ1bmN0aW9uKHRhYikge1xuICBmb3IodmFyIGk9MDsgaTwyNTY7ICsraSkge1xuICAgIHZhciB2ID0gaSwgciA9IGksIHMgPSA3O1xuICAgIGZvciAodiA+Pj49IDE7IHY7IHYgPj4+PSAxKSB7XG4gICAgICByIDw8PSAxO1xuICAgICAgciB8PSB2ICYgMTtcbiAgICAgIC0tcztcbiAgICB9XG4gICAgdGFiW2ldID0gKHIgPDwgcykgJiAweGZmO1xuICB9XG59KShSRVZFUlNFX1RBQkxFKTtcblxuLy9SZXZlcnNlIGJpdHMgaW4gYSAzMiBiaXQgd29yZFxuZXhwb3J0cy5yZXZlcnNlID0gZnVuY3Rpb24odikge1xuICByZXR1cm4gIChSRVZFUlNFX1RBQkxFWyB2ICAgICAgICAgJiAweGZmXSA8PCAyNCkgfFxuICAgICAgICAgIChSRVZFUlNFX1RBQkxFWyh2ID4+PiA4KSAgJiAweGZmXSA8PCAxNikgfFxuICAgICAgICAgIChSRVZFUlNFX1RBQkxFWyh2ID4+PiAxNikgJiAweGZmXSA8PCA4KSAgfFxuICAgICAgICAgICBSRVZFUlNFX1RBQkxFWyh2ID4+PiAyNCkgJiAweGZmXTtcbn1cblxuLy9JbnRlcmxlYXZlIGJpdHMgb2YgMiBjb29yZGluYXRlcyB3aXRoIDE2IGJpdHMuICBVc2VmdWwgZm9yIGZhc3QgcXVhZHRyZWUgY29kZXNcbmV4cG9ydHMuaW50ZXJsZWF2ZTIgPSBmdW5jdGlvbih4LCB5KSB7XG4gIHggJj0gMHhGRkZGO1xuICB4ID0gKHggfCAoeCA8PCA4KSkgJiAweDAwRkYwMEZGO1xuICB4ID0gKHggfCAoeCA8PCA0KSkgJiAweDBGMEYwRjBGO1xuICB4ID0gKHggfCAoeCA8PCAyKSkgJiAweDMzMzMzMzMzO1xuICB4ID0gKHggfCAoeCA8PCAxKSkgJiAweDU1NTU1NTU1O1xuXG4gIHkgJj0gMHhGRkZGO1xuICB5ID0gKHkgfCAoeSA8PCA4KSkgJiAweDAwRkYwMEZGO1xuICB5ID0gKHkgfCAoeSA8PCA0KSkgJiAweDBGMEYwRjBGO1xuICB5ID0gKHkgfCAoeSA8PCAyKSkgJiAweDMzMzMzMzMzO1xuICB5ID0gKHkgfCAoeSA8PCAxKSkgJiAweDU1NTU1NTU1O1xuXG4gIHJldHVybiB4IHwgKHkgPDwgMSk7XG59XG5cbi8vRXh0cmFjdHMgdGhlIG50aCBpbnRlcmxlYXZlZCBjb21wb25lbnRcbmV4cG9ydHMuZGVpbnRlcmxlYXZlMiA9IGZ1bmN0aW9uKHYsIG4pIHtcbiAgdiA9ICh2ID4+PiBuKSAmIDB4NTU1NTU1NTU7XG4gIHYgPSAodiB8ICh2ID4+PiAxKSkgICYgMHgzMzMzMzMzMztcbiAgdiA9ICh2IHwgKHYgPj4+IDIpKSAgJiAweDBGMEYwRjBGO1xuICB2ID0gKHYgfCAodiA+Pj4gNCkpICAmIDB4MDBGRjAwRkY7XG4gIHYgPSAodiB8ICh2ID4+PiAxNikpICYgMHgwMDBGRkZGO1xuICByZXR1cm4gKHYgPDwgMTYpID4+IDE2O1xufVxuXG5cbi8vSW50ZXJsZWF2ZSBiaXRzIG9mIDMgY29vcmRpbmF0ZXMsIGVhY2ggd2l0aCAxMCBiaXRzLiAgVXNlZnVsIGZvciBmYXN0IG9jdHJlZSBjb2Rlc1xuZXhwb3J0cy5pbnRlcmxlYXZlMyA9IGZ1bmN0aW9uKHgsIHksIHopIHtcbiAgeCAmPSAweDNGRjtcbiAgeCAgPSAoeCB8ICh4PDwxNikpICYgNDI3ODE5MDMzNTtcbiAgeCAgPSAoeCB8ICh4PDw4KSkgICYgMjUxNzE5Njk1O1xuICB4ICA9ICh4IHwgKHg8PDQpKSAgJiAzMjcyMzU2MDM1O1xuICB4ICA9ICh4IHwgKHg8PDIpKSAgJiAxMjI3MTMzNTEzO1xuXG4gIHkgJj0gMHgzRkY7XG4gIHkgID0gKHkgfCAoeTw8MTYpKSAmIDQyNzgxOTAzMzU7XG4gIHkgID0gKHkgfCAoeTw8OCkpICAmIDI1MTcxOTY5NTtcbiAgeSAgPSAoeSB8ICh5PDw0KSkgICYgMzI3MjM1NjAzNTtcbiAgeSAgPSAoeSB8ICh5PDwyKSkgICYgMTIyNzEzMzUxMztcbiAgeCB8PSAoeSA8PCAxKTtcbiAgXG4gIHogJj0gMHgzRkY7XG4gIHogID0gKHogfCAoejw8MTYpKSAmIDQyNzgxOTAzMzU7XG4gIHogID0gKHogfCAoejw8OCkpICAmIDI1MTcxOTY5NTtcbiAgeiAgPSAoeiB8ICh6PDw0KSkgICYgMzI3MjM1NjAzNTtcbiAgeiAgPSAoeiB8ICh6PDwyKSkgICYgMTIyNzEzMzUxMztcbiAgXG4gIHJldHVybiB4IHwgKHogPDwgMik7XG59XG5cbi8vRXh0cmFjdHMgbnRoIGludGVybGVhdmVkIGNvbXBvbmVudCBvZiBhIDMtdHVwbGVcbmV4cG9ydHMuZGVpbnRlcmxlYXZlMyA9IGZ1bmN0aW9uKHYsIG4pIHtcbiAgdiA9ICh2ID4+PiBuKSAgICAgICAmIDEyMjcxMzM1MTM7XG4gIHYgPSAodiB8ICh2Pj4+MikpICAgJiAzMjcyMzU2MDM1O1xuICB2ID0gKHYgfCAodj4+PjQpKSAgICYgMjUxNzE5Njk1O1xuICB2ID0gKHYgfCAodj4+PjgpKSAgICYgNDI3ODE5MDMzNTtcbiAgdiA9ICh2IHwgKHY+Pj4xNikpICAmIDB4M0ZGO1xuICByZXR1cm4gKHY8PDIyKT4+MjI7XG59XG5cbi8vQ29tcHV0ZXMgbmV4dCBjb21iaW5hdGlvbiBpbiBjb2xleGljb2dyYXBoaWMgb3JkZXIgKHRoaXMgaXMgbWlzdGFrZW5seSBjYWxsZWQgbmV4dFBlcm11dGF0aW9uIG9uIHRoZSBiaXQgdHdpZGRsaW5nIGhhY2tzIHBhZ2UpXG5leHBvcnRzLm5leHRDb21iaW5hdGlvbiA9IGZ1bmN0aW9uKHYpIHtcbiAgdmFyIHQgPSB2IHwgKHYgLSAxKTtcbiAgcmV0dXJuICh0ICsgMSkgfCAoKCh+dCAmIC1+dCkgLSAxKSA+Pj4gKGNvdW50VHJhaWxpbmdaZXJvcyh2KSArIDEpKTtcbn1cblxuIiwiXCJ1c2Ugc3RyaWN0XCI7IFwidXNlIHJlc3RyaWN0XCI7XG5cbm1vZHVsZS5leHBvcnRzID0gVW5pb25GaW5kO1xuXG5mdW5jdGlvbiBVbmlvbkZpbmQoY291bnQpIHtcbiAgdGhpcy5yb290cyA9IG5ldyBBcnJheShjb3VudCk7XG4gIHRoaXMucmFua3MgPSBuZXcgQXJyYXkoY291bnQpO1xuICBcbiAgZm9yKHZhciBpPTA7IGk8Y291bnQ7ICsraSkge1xuICAgIHRoaXMucm9vdHNbaV0gPSBpO1xuICAgIHRoaXMucmFua3NbaV0gPSAwO1xuICB9XG59XG5cbnZhciBwcm90byA9IFVuaW9uRmluZC5wcm90b3R5cGVcblxuT2JqZWN0LmRlZmluZVByb3BlcnR5KHByb3RvLCBcImxlbmd0aFwiLCB7XG4gIFwiZ2V0XCI6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnJvb3RzLmxlbmd0aFxuICB9XG59KVxuXG5wcm90by5tYWtlU2V0ID0gZnVuY3Rpb24oKSB7XG4gIHZhciBuID0gdGhpcy5yb290cy5sZW5ndGg7XG4gIHRoaXMucm9vdHMucHVzaChuKTtcbiAgdGhpcy5yYW5rcy5wdXNoKDApO1xuICByZXR1cm4gbjtcbn1cblxucHJvdG8uZmluZCA9IGZ1bmN0aW9uKHgpIHtcbiAgdmFyIHJvb3RzID0gdGhpcy5yb290cztcbiAgd2hpbGUocm9vdHNbeF0gIT09IHgpIHtcbiAgICB2YXIgeSA9IHJvb3RzW3hdO1xuICAgIHJvb3RzW3hdID0gcm9vdHNbeV07XG4gICAgeCA9IHk7XG4gIH1cbiAgcmV0dXJuIHg7XG59XG5cbnByb3RvLmxpbmsgPSBmdW5jdGlvbih4LCB5KSB7XG4gIHZhciB4ciA9IHRoaXMuZmluZCh4KVxuICAgICwgeXIgPSB0aGlzLmZpbmQoeSk7XG4gIGlmKHhyID09PSB5cikge1xuICAgIHJldHVybjtcbiAgfVxuICB2YXIgcmFua3MgPSB0aGlzLnJhbmtzXG4gICAgLCByb290cyA9IHRoaXMucm9vdHNcbiAgICAsIHhkICAgID0gcmFua3NbeHJdXG4gICAgLCB5ZCAgICA9IHJhbmtzW3lyXTtcbiAgaWYoeGQgPCB5ZCkge1xuICAgIHJvb3RzW3hyXSA9IHlyO1xuICB9IGVsc2UgaWYoeWQgPCB4ZCkge1xuICAgIHJvb3RzW3lyXSA9IHhyO1xuICB9IGVsc2Uge1xuICAgIHJvb3RzW3lyXSA9IHhyO1xuICAgICsrcmFua3NbeHJdO1xuICB9XG59IiwiXCJ1c2Ugc3RyaWN0XCI7IFwidXNlIHJlc3RyaWN0XCI7XG5cbnZhciBiaXRzICAgICAgPSByZXF1aXJlKFwiYml0LXR3aWRkbGVcIilcbiAgLCBVbmlvbkZpbmQgPSByZXF1aXJlKFwidW5pb24tZmluZFwiKVxuXG4vL1JldHVybnMgdGhlIGRpbWVuc2lvbiBvZiBhIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gZGltZW5zaW9uKGNlbGxzKSB7XG4gIHZhciBkID0gMFxuICAgICwgbWF4ID0gTWF0aC5tYXhcbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgZCA9IG1heChkLCBjZWxsc1tpXS5sZW5ndGgpXG4gIH1cbiAgcmV0dXJuIGQtMVxufVxuZXhwb3J0cy5kaW1lbnNpb24gPSBkaW1lbnNpb25cblxuLy9Db3VudHMgdGhlIG51bWJlciBvZiB2ZXJ0aWNlcyBpbiBmYWNlc1xuZnVuY3Rpb24gY291bnRWZXJ0aWNlcyhjZWxscykge1xuICB2YXIgdmMgPSAtMVxuICAgICwgbWF4ID0gTWF0aC5tYXhcbiAgZm9yKHZhciBpPTAsIGlsPWNlbGxzLmxlbmd0aDsgaTxpbDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wLCBqbD1jLmxlbmd0aDsgajxqbDsgKytqKSB7XG4gICAgICB2YyA9IG1heCh2YywgY1tqXSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHZjKzFcbn1cbmV4cG9ydHMuY291bnRWZXJ0aWNlcyA9IGNvdW50VmVydGljZXNcblxuLy9SZXR1cm5zIGEgZGVlcCBjb3B5IG9mIGNlbGxzXG5mdW5jdGlvbiBjbG9uZUNlbGxzKGNlbGxzKSB7XG4gIHZhciBuY2VsbHMgPSBuZXcgQXJyYXkoY2VsbHMubGVuZ3RoKVxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICBuY2VsbHNbaV0gPSBjZWxsc1tpXS5zbGljZSgwKVxuICB9XG4gIHJldHVybiBuY2VsbHNcbn1cbmV4cG9ydHMuY2xvbmVDZWxscyA9IGNsb25lQ2VsbHNcblxuLy9SYW5rcyBhIHBhaXIgb2YgY2VsbHMgdXAgdG8gcGVybXV0YXRpb25cbmZ1bmN0aW9uIGNvbXBhcmVDZWxscyhhLCBiKSB7XG4gIHZhciBuID0gYS5sZW5ndGhcbiAgICAsIHQgPSBhLmxlbmd0aCAtIGIubGVuZ3RoXG4gICAgLCBtaW4gPSBNYXRoLm1pblxuICBpZih0KSB7XG4gICAgcmV0dXJuIHRcbiAgfVxuICBzd2l0Y2gobikge1xuICAgIGNhc2UgMDpcbiAgICAgIHJldHVybiAwO1xuICAgIGNhc2UgMTpcbiAgICAgIHJldHVybiBhWzBdIC0gYlswXTtcbiAgICBjYXNlIDI6XG4gICAgICB2YXIgZCA9IGFbMF0rYVsxXS1iWzBdLWJbMV1cbiAgICAgIGlmKGQpIHtcbiAgICAgICAgcmV0dXJuIGRcbiAgICAgIH1cbiAgICAgIHJldHVybiBtaW4oYVswXSxhWzFdKSAtIG1pbihiWzBdLGJbMV0pXG4gICAgY2FzZSAzOlxuICAgICAgdmFyIGwxID0gYVswXSthWzFdXG4gICAgICAgICwgbTEgPSBiWzBdK2JbMV1cbiAgICAgIGQgPSBsMSthWzJdIC0gKG0xK2JbMl0pXG4gICAgICBpZihkKSB7XG4gICAgICAgIHJldHVybiBkXG4gICAgICB9XG4gICAgICB2YXIgbDAgPSBtaW4oYVswXSwgYVsxXSlcbiAgICAgICAgLCBtMCA9IG1pbihiWzBdLCBiWzFdKVxuICAgICAgICAsIGQgID0gbWluKGwwLCBhWzJdKSAtIG1pbihtMCwgYlsyXSlcbiAgICAgIGlmKGQpIHtcbiAgICAgICAgcmV0dXJuIGRcbiAgICAgIH1cbiAgICAgIHJldHVybiBtaW4obDArYVsyXSwgbDEpIC0gbWluKG0wK2JbMl0sIG0xKVxuICAgIFxuICAgIC8vVE9ETzogTWF5YmUgb3B0aW1pemUgbj00IGFzIHdlbGw/XG4gICAgXG4gICAgZGVmYXVsdDpcbiAgICAgIHZhciBhcyA9IGEuc2xpY2UoMClcbiAgICAgIGFzLnNvcnQoKVxuICAgICAgdmFyIGJzID0gYi5zbGljZSgwKVxuICAgICAgYnMuc29ydCgpXG4gICAgICBmb3IodmFyIGk9MDsgaTxuOyArK2kpIHtcbiAgICAgICAgdCA9IGFzW2ldIC0gYnNbaV1cbiAgICAgICAgaWYodCkge1xuICAgICAgICAgIHJldHVybiB0XG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJldHVybiAwXG4gIH1cbn1cbmV4cG9ydHMuY29tcGFyZUNlbGxzID0gY29tcGFyZUNlbGxzXG5cbmZ1bmN0aW9uIGNvbXBhcmVaaXBwZWQoYSwgYikge1xuICByZXR1cm4gY29tcGFyZUNlbGxzKGFbMF0sIGJbMF0pXG59XG5cbi8vUHV0cyBhIGNlbGwgY29tcGxleCBpbnRvIG5vcm1hbCBvcmRlciBmb3IgdGhlIHB1cnBvc2VzIG9mIGZpbmRDZWxsIHF1ZXJpZXNcbmZ1bmN0aW9uIG5vcm1hbGl6ZShjZWxscywgYXR0cikge1xuICBpZihhdHRyKSB7XG4gICAgdmFyIGxlbiA9IGNlbGxzLmxlbmd0aFxuICAgIHZhciB6aXBwZWQgPSBuZXcgQXJyYXkobGVuKVxuICAgIGZvcih2YXIgaT0wOyBpPGxlbjsgKytpKSB7XG4gICAgICB6aXBwZWRbaV0gPSBbY2VsbHNbaV0sIGF0dHJbaV1dXG4gICAgfVxuICAgIHppcHBlZC5zb3J0KGNvbXBhcmVaaXBwZWQpXG4gICAgZm9yKHZhciBpPTA7IGk8bGVuOyArK2kpIHtcbiAgICAgIGNlbGxzW2ldID0gemlwcGVkW2ldWzBdXG4gICAgICBhdHRyW2ldID0gemlwcGVkW2ldWzFdXG4gICAgfVxuICAgIHJldHVybiBjZWxsc1xuICB9IGVsc2Uge1xuICAgIGNlbGxzLnNvcnQoY29tcGFyZUNlbGxzKVxuICAgIHJldHVybiBjZWxsc1xuICB9XG59XG5leHBvcnRzLm5vcm1hbGl6ZSA9IG5vcm1hbGl6ZVxuXG4vL1JlbW92ZXMgYWxsIGR1cGxpY2F0ZSBjZWxscyBpbiB0aGUgY29tcGxleFxuZnVuY3Rpb24gdW5pcXVlKGNlbGxzKSB7XG4gIGlmKGNlbGxzLmxlbmd0aCA9PT0gMCkge1xuICAgIHJldHVybiBbXVxuICB9XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBjZWxscy5sZW5ndGhcbiAgZm9yKHZhciBpPTE7IGk8bGVuOyArK2kpIHtcbiAgICB2YXIgYSA9IGNlbGxzW2ldXG4gICAgaWYoY29tcGFyZUNlbGxzKGEsIGNlbGxzW2ktMV0pKSB7XG4gICAgICBpZihpID09PSBwdHIpIHtcbiAgICAgICAgcHRyKytcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGNlbGxzW3B0cisrXSA9IGFcbiAgICB9XG4gIH1cbiAgY2VsbHMubGVuZ3RoID0gcHRyXG4gIHJldHVybiBjZWxsc1xufVxuZXhwb3J0cy51bmlxdWUgPSB1bmlxdWU7XG5cbi8vRmluZHMgYSBjZWxsIGluIGEgbm9ybWFsaXplZCBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGZpbmRDZWxsKGNlbGxzLCBjKSB7XG4gIHZhciBsbyA9IDBcbiAgICAsIGhpID0gY2VsbHMubGVuZ3RoLTFcbiAgICAsIHIgID0gLTFcbiAgd2hpbGUgKGxvIDw9IGhpKSB7XG4gICAgdmFyIG1pZCA9IChsbyArIGhpKSA+PiAxXG4gICAgICAsIHMgICA9IGNvbXBhcmVDZWxscyhjZWxsc1ttaWRdLCBjKVxuICAgIGlmKHMgPD0gMCkge1xuICAgICAgaWYocyA9PT0gMCkge1xuICAgICAgICByID0gbWlkXG4gICAgICB9XG4gICAgICBsbyA9IG1pZCArIDFcbiAgICB9IGVsc2UgaWYocyA+IDApIHtcbiAgICAgIGhpID0gbWlkIC0gMVxuICAgIH1cbiAgfVxuICByZXR1cm4gclxufVxuZXhwb3J0cy5maW5kQ2VsbCA9IGZpbmRDZWxsO1xuXG4vL0J1aWxkcyBhbiBpbmRleCBmb3IgYW4gbi1jZWxsLiAgVGhpcyBpcyBtb3JlIGdlbmVyYWwgdGhhbiBkdWFsLCBidXQgbGVzcyBlZmZpY2llbnRcbmZ1bmN0aW9uIGluY2lkZW5jZShmcm9tX2NlbGxzLCB0b19jZWxscykge1xuICB2YXIgaW5kZXggPSBuZXcgQXJyYXkoZnJvbV9jZWxscy5sZW5ndGgpXG4gIGZvcih2YXIgaT0wLCBpbD1pbmRleC5sZW5ndGg7IGk8aWw7ICsraSkge1xuICAgIGluZGV4W2ldID0gW11cbiAgfVxuICB2YXIgYiA9IFtdXG4gIGZvcih2YXIgaT0wLCBuPXRvX2NlbGxzLmxlbmd0aDsgaTxuOyArK2kpIHtcbiAgICB2YXIgYyA9IHRvX2NlbGxzW2ldXG4gICAgdmFyIGNsID0gYy5sZW5ndGhcbiAgICBmb3IodmFyIGs9MSwga249KDE8PGNsKTsgazxrbjsgKytrKSB7XG4gICAgICBiLmxlbmd0aCA9IGJpdHMucG9wQ291bnQoaylcbiAgICAgIHZhciBsID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8Y2w7ICsraikge1xuICAgICAgICBpZihrICYgKDE8PGopKSB7XG4gICAgICAgICAgYltsKytdID0gY1tqXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICB2YXIgaWR4PWZpbmRDZWxsKGZyb21fY2VsbHMsIGIpXG4gICAgICBpZihpZHggPCAwKSB7XG4gICAgICAgIGNvbnRpbnVlXG4gICAgICB9XG4gICAgICB3aGlsZSh0cnVlKSB7XG4gICAgICAgIGluZGV4W2lkeCsrXS5wdXNoKGkpXG4gICAgICAgIGlmKGlkeCA+PSBmcm9tX2NlbGxzLmxlbmd0aCB8fCBjb21wYXJlQ2VsbHMoZnJvbV9jZWxsc1tpZHhdLCBiKSAhPT0gMCkge1xuICAgICAgICAgIGJyZWFrXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gIH1cbiAgcmV0dXJuIGluZGV4XG59XG5leHBvcnRzLmluY2lkZW5jZSA9IGluY2lkZW5jZVxuXG4vL0NvbXB1dGVzIHRoZSBkdWFsIG9mIHRoZSBtZXNoLiAgVGhpcyBpcyBiYXNpY2FsbHkgYW4gb3B0aW1pemVkIHZlcnNpb24gb2YgYnVpbGRJbmRleCBmb3IgdGhlIHNpdHVhdGlvbiB3aGVyZSBmcm9tX2NlbGxzIGlzIGp1c3QgdGhlIGxpc3Qgb2YgdmVydGljZXNcbmZ1bmN0aW9uIGR1YWwoY2VsbHMsIHZlcnRleF9jb3VudCkge1xuICBpZighdmVydGV4X2NvdW50KSB7XG4gICAgcmV0dXJuIGluY2lkZW5jZSh1bmlxdWUoc2tlbGV0b24oY2VsbHMsIDApKSwgY2VsbHMsIDApXG4gIH1cbiAgdmFyIHJlcyA9IG5ldyBBcnJheSh2ZXJ0ZXhfY291bnQpXG4gIGZvcih2YXIgaT0wOyBpPHZlcnRleF9jb3VudDsgKytpKSB7XG4gICAgcmVzW2ldID0gW11cbiAgfVxuICBmb3IodmFyIGk9MCwgbGVuPWNlbGxzLmxlbmd0aDsgaTxsZW47ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGo9MCwgY2w9Yy5sZW5ndGg7IGo8Y2w7ICsraikge1xuICAgICAgcmVzW2Nbal1dLnB1c2goaSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIHJlc1xufVxuZXhwb3J0cy5kdWFsID0gZHVhbFxuXG4vL0VudW1lcmF0ZXMgYWxsIGNlbGxzIGluIHRoZSBjb21wbGV4XG5mdW5jdGlvbiBleHBsb2RlKGNlbGxzKSB7XG4gIHZhciByZXN1bHQgPSBbXVxuICBmb3IodmFyIGk9MCwgaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgICAsIGNsID0gYy5sZW5ndGh8MFxuICAgIGZvcih2YXIgaj0xLCBqbD0oMTw8Y2wpOyBqPGpsOyArK2opIHtcbiAgICAgIHZhciBiID0gW11cbiAgICAgIGZvcih2YXIgaz0wOyBrPGNsOyArK2spIHtcbiAgICAgICAgaWYoKGogPj4+IGspICYgMSkge1xuICAgICAgICAgIGIucHVzaChjW2tdKVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXN1bHQucHVzaChiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gbm9ybWFsaXplKHJlc3VsdClcbn1cbmV4cG9ydHMuZXhwbG9kZSA9IGV4cGxvZGVcblxuLy9FbnVtZXJhdGVzIGFsbCBvZiB0aGUgbi1jZWxscyBvZiBhIGNlbGwgY29tcGxleFxuZnVuY3Rpb24gc2tlbGV0b24oY2VsbHMsIG4pIHtcbiAgaWYobiA8IDApIHtcbiAgICByZXR1cm4gW11cbiAgfVxuICB2YXIgcmVzdWx0ID0gW11cbiAgICAsIGswICAgICA9ICgxPDwobisxKSktMVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBjID0gY2VsbHNbaV1cbiAgICBmb3IodmFyIGs9azA7IGs8KDE8PGMubGVuZ3RoKTsgaz1iaXRzLm5leHRDb21iaW5hdGlvbihrKSkge1xuICAgICAgdmFyIGIgPSBuZXcgQXJyYXkobisxKVxuICAgICAgICAsIGwgPSAwXG4gICAgICBmb3IodmFyIGo9MDsgajxjLmxlbmd0aDsgKytqKSB7XG4gICAgICAgIGlmKGsgJiAoMTw8aikpIHtcbiAgICAgICAgICBiW2wrK10gPSBjW2pdXG4gICAgICAgIH1cbiAgICAgIH1cbiAgICAgIHJlc3VsdC5wdXNoKGIpXG4gICAgfVxuICB9XG4gIHJldHVybiBub3JtYWxpemUocmVzdWx0KVxufVxuZXhwb3J0cy5za2VsZXRvbiA9IHNrZWxldG9uO1xuXG4vL0NvbXB1dGVzIHRoZSBib3VuZGFyeSBvZiBhbGwgY2VsbHMsIGRvZXMgbm90IHJlbW92ZSBkdXBsaWNhdGVzXG5mdW5jdGlvbiBib3VuZGFyeShjZWxscykge1xuICB2YXIgcmVzID0gW11cbiAgZm9yKHZhciBpPTAsaWw9Y2VsbHMubGVuZ3RoOyBpPGlsOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTAsY2w9Yy5sZW5ndGg7IGo8Y2w7ICsraikge1xuICAgICAgdmFyIGIgPSBuZXcgQXJyYXkoYy5sZW5ndGgtMSlcbiAgICAgIGZvcih2YXIgaz0wLCBsPTA7IGs8Y2w7ICsraykge1xuICAgICAgICBpZihrICE9PSBqKSB7XG4gICAgICAgICAgYltsKytdID0gY1trXVxuICAgICAgICB9XG4gICAgICB9XG4gICAgICByZXMucHVzaChiKVxuICAgIH1cbiAgfVxuICByZXR1cm4gbm9ybWFsaXplKHJlcylcbn1cbmV4cG9ydHMuYm91bmRhcnkgPSBib3VuZGFyeTtcblxuLy9Db21wdXRlcyBjb25uZWN0ZWQgY29tcG9uZW50cyBmb3IgYSBkZW5zZSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGNvbm5lY3RlZENvbXBvbmVudHNfZGVuc2UoY2VsbHMsIHZlcnRleF9jb3VudCkge1xuICB2YXIgbGFiZWxzID0gbmV3IFVuaW9uRmluZCh2ZXJ0ZXhfY291bnQpXG4gIGZvcih2YXIgaT0wOyBpPGNlbGxzLmxlbmd0aDsgKytpKSB7XG4gICAgdmFyIGMgPSBjZWxsc1tpXVxuICAgIGZvcih2YXIgaj0wOyBqPGMubGVuZ3RoOyArK2opIHtcbiAgICAgIGZvcih2YXIgaz1qKzE7IGs8Yy5sZW5ndGg7ICsraykge1xuICAgICAgICBsYWJlbHMubGluayhjW2pdLCBjW2tdKVxuICAgICAgfVxuICAgIH1cbiAgfVxuICB2YXIgY29tcG9uZW50cyA9IFtdXG4gICAgLCBjb21wb25lbnRfbGFiZWxzID0gbGFiZWxzLnJhbmtzXG4gIGZvcih2YXIgaT0wOyBpPGNvbXBvbmVudF9sYWJlbHMubGVuZ3RoOyArK2kpIHtcbiAgICBjb21wb25lbnRfbGFiZWxzW2ldID0gLTFcbiAgfVxuICBmb3IodmFyIGk9MDsgaTxjZWxscy5sZW5ndGg7ICsraSkge1xuICAgIHZhciBsID0gbGFiZWxzLmZpbmQoY2VsbHNbaV1bMF0pXG4gICAgaWYoY29tcG9uZW50X2xhYmVsc1tsXSA8IDApIHtcbiAgICAgIGNvbXBvbmVudF9sYWJlbHNbbF0gPSBjb21wb25lbnRzLmxlbmd0aFxuICAgICAgY29tcG9uZW50cy5wdXNoKFtjZWxsc1tpXS5zbGljZSgwKV0pXG4gICAgfSBlbHNlIHtcbiAgICAgIGNvbXBvbmVudHNbY29tcG9uZW50X2xhYmVsc1tsXV0ucHVzaChjZWxsc1tpXS5zbGljZSgwKSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIGNvbXBvbmVudHNcbn1cblxuLy9Db21wdXRlcyBjb25uZWN0ZWQgY29tcG9uZW50cyBmb3IgYSBzcGFyc2UgZ3JhcGhcbmZ1bmN0aW9uIGNvbm5lY3RlZENvbXBvbmVudHNfc3BhcnNlKGNlbGxzKSB7XG4gIHZhciB2ZXJ0aWNlcyAgPSB1bmlxdWUobm9ybWFsaXplKHNrZWxldG9uKGNlbGxzLCAwKSkpXG4gICAgLCBsYWJlbHMgICAgPSBuZXcgVW5pb25GaW5kKHZlcnRpY2VzLmxlbmd0aClcbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgYyA9IGNlbGxzW2ldXG4gICAgZm9yKHZhciBqPTA7IGo8Yy5sZW5ndGg7ICsraikge1xuICAgICAgdmFyIHZqID0gZmluZENlbGwodmVydGljZXMsIFtjW2pdXSlcbiAgICAgIGZvcih2YXIgaz1qKzE7IGs8Yy5sZW5ndGg7ICsraykge1xuICAgICAgICBsYWJlbHMubGluayh2aiwgZmluZENlbGwodmVydGljZXMsIFtjW2tdXSkpXG4gICAgICB9XG4gICAgfVxuICB9XG4gIHZhciBjb21wb25lbnRzICAgICAgICA9IFtdXG4gICAgLCBjb21wb25lbnRfbGFiZWxzICA9IGxhYmVscy5yYW5rc1xuICBmb3IodmFyIGk9MDsgaTxjb21wb25lbnRfbGFiZWxzLmxlbmd0aDsgKytpKSB7XG4gICAgY29tcG9uZW50X2xhYmVsc1tpXSA9IC0xXG4gIH1cbiAgZm9yKHZhciBpPTA7IGk8Y2VsbHMubGVuZ3RoOyArK2kpIHtcbiAgICB2YXIgbCA9IGxhYmVscy5maW5kKGZpbmRDZWxsKHZlcnRpY2VzLCBbY2VsbHNbaV1bMF1dKSk7XG4gICAgaWYoY29tcG9uZW50X2xhYmVsc1tsXSA8IDApIHtcbiAgICAgIGNvbXBvbmVudF9sYWJlbHNbbF0gPSBjb21wb25lbnRzLmxlbmd0aFxuICAgICAgY29tcG9uZW50cy5wdXNoKFtjZWxsc1tpXS5zbGljZSgwKV0pXG4gICAgfSBlbHNlIHtcbiAgICAgIGNvbXBvbmVudHNbY29tcG9uZW50X2xhYmVsc1tsXV0ucHVzaChjZWxsc1tpXS5zbGljZSgwKSlcbiAgICB9XG4gIH1cbiAgcmV0dXJuIGNvbXBvbmVudHNcbn1cblxuLy9Db21wdXRlcyBjb25uZWN0ZWQgY29tcG9uZW50cyBmb3IgYSBjZWxsIGNvbXBsZXhcbmZ1bmN0aW9uIGNvbm5lY3RlZENvbXBvbmVudHMoY2VsbHMsIHZlcnRleF9jb3VudCkge1xuICBpZih2ZXJ0ZXhfY291bnQpIHtcbiAgICByZXR1cm4gY29ubmVjdGVkQ29tcG9uZW50c19kZW5zZShjZWxscywgdmVydGV4X2NvdW50KVxuICB9XG4gIHJldHVybiBjb25uZWN0ZWRDb21wb25lbnRzX3NwYXJzZShjZWxscylcbn1cbmV4cG9ydHMuY29ubmVjdGVkQ29tcG9uZW50cyA9IGNvbm5lY3RlZENvbXBvbmVudHNcbiIsIlwidXNlIHN0cmljdFwiXG5cbmZ1bmN0aW9uIHVuaXF1ZV9wcmVkKGxpc3QsIGNvbXBhcmUpIHtcbiAgdmFyIHB0ciA9IDFcbiAgICAsIGxlbiA9IGxpc3QubGVuZ3RoXG4gICAgLCBhPWxpc3RbMF0sIGI9bGlzdFswXVxuICBmb3IodmFyIGk9MTsgaTxsZW47ICsraSkge1xuICAgIGIgPSBhXG4gICAgYSA9IGxpc3RbaV1cbiAgICBpZihjb21wYXJlKGEsIGIpKSB7XG4gICAgICBpZihpID09PSBwdHIpIHtcbiAgICAgICAgcHRyKytcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGxpc3RbcHRyKytdID0gYVxuICAgIH1cbiAgfVxuICBsaXN0Lmxlbmd0aCA9IHB0clxuICByZXR1cm4gbGlzdFxufVxuXG5mdW5jdGlvbiB1bmlxdWVfZXEobGlzdCkge1xuICB2YXIgcHRyID0gMVxuICAgICwgbGVuID0gbGlzdC5sZW5ndGhcbiAgICAsIGE9bGlzdFswXSwgYiA9IGxpc3RbMF1cbiAgZm9yKHZhciBpPTE7IGk8bGVuOyArK2ksIGI9YSkge1xuICAgIGIgPSBhXG4gICAgYSA9IGxpc3RbaV1cbiAgICBpZihhICE9PSBiKSB7XG4gICAgICBpZihpID09PSBwdHIpIHtcbiAgICAgICAgcHRyKytcbiAgICAgICAgY29udGludWVcbiAgICAgIH1cbiAgICAgIGxpc3RbcHRyKytdID0gYVxuICAgIH1cbiAgfVxuICBsaXN0Lmxlbmd0aCA9IHB0clxuICByZXR1cm4gbGlzdFxufVxuXG5mdW5jdGlvbiB1bmlxdWUobGlzdCwgY29tcGFyZSwgc29ydGVkKSB7XG4gIGlmKGxpc3QubGVuZ3RoID09PSAwKSB7XG4gICAgcmV0dXJuIGxpc3RcbiAgfVxuICBpZihjb21wYXJlKSB7XG4gICAgaWYoIXNvcnRlZCkge1xuICAgICAgbGlzdC5zb3J0KGNvbXBhcmUpXG4gICAgfVxuICAgIHJldHVybiB1bmlxdWVfcHJlZChsaXN0LCBjb21wYXJlKVxuICB9XG4gIGlmKCFzb3J0ZWQpIHtcbiAgICBsaXN0LnNvcnQoKVxuICB9XG4gIHJldHVybiB1bmlxdWVfZXEobGlzdClcbn1cblxubW9kdWxlLmV4cG9ydHMgPSB1bmlxdWVcbiIsIlwidXNlIHN0cmljdFwiXG5cbnZhciBjaCA9IHJlcXVpcmUoXCJpbmNyZW1lbnRhbC1jb252ZXgtaHVsbFwiKVxudmFyIHVuaXEgPSByZXF1aXJlKFwidW5pcVwiKVxuXG5tb2R1bGUuZXhwb3J0cyA9IHRyaWFuZ3VsYXRlXG5cbmZ1bmN0aW9uIExpZnRlZFBvaW50KHAsIGkpIHtcbiAgdGhpcy5wb2ludCA9IHBcbiAgdGhpcy5pbmRleCA9IGlcbn1cblxuZnVuY3Rpb24gY29tcGFyZUxpZnRlZChhLCBiKSB7XG4gIHZhciBhcCA9IGEucG9pbnRcbiAgdmFyIGJwID0gYi5wb2ludFxuICB2YXIgZCA9IGFwLmxlbmd0aFxuICBmb3IodmFyIGk9MDsgaTxkOyArK2kpIHtcbiAgICB2YXIgcyA9IGJwW2ldIC0gYXBbaV1cbiAgICBpZihzKSB7XG4gICAgICByZXR1cm4gc1xuICAgIH1cbiAgfVxuICByZXR1cm4gMFxufVxuXG5mdW5jdGlvbiB0cmlhbmd1bGF0ZTFEKG4sIHBvaW50cywgaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICBpZihuID09PSAxKSB7XG4gICAgaWYoaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICAgICAgcmV0dXJuIFsgWy0xLCAwXSBdXG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBbXVxuICAgIH1cbiAgfVxuICB2YXIgbGlmdGVkID0gcG9pbnRzLm1hcChmdW5jdGlvbihwLCBpKSB7XG4gICAgcmV0dXJuIFsgcFswXSwgaSBdXG4gIH0pXG4gIGxpZnRlZC5zb3J0KGZ1bmN0aW9uKGEsYikge1xuICAgIHJldHVybiBhWzBdIC0gYlswXVxuICB9KVxuICB2YXIgY2VsbHMgPSBuZXcgQXJyYXkobiAtIDEpXG4gIGZvcih2YXIgaT0xOyBpPG47ICsraSkge1xuICAgIHZhciBhID0gbGlmdGVkW2ktMV1cbiAgICB2YXIgYiA9IGxpZnRlZFtpXVxuICAgIGNlbGxzW2ktMV0gPSBbIGFbMV0sIGJbMV0gXVxuICB9XG4gIGlmKGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgICBjZWxscy5wdXNoKFxuICAgICAgWyAtMSwgY2VsbHNbMF1bMV0sIF0sXG4gICAgICBbIGNlbGxzW24tMV1bMV0sIC0xIF0pXG4gIH1cbiAgcmV0dXJuIGNlbGxzXG59XG5cbmZ1bmN0aW9uIHRyaWFuZ3VsYXRlKHBvaW50cywgaW5jbHVkZVBvaW50QXRJbmZpbml0eSkge1xuICB2YXIgbiA9IHBvaW50cy5sZW5ndGhcbiAgaWYobiA9PT0gMCkge1xuICAgIHJldHVybiBbXVxuICB9XG4gIFxuICB2YXIgZCA9IHBvaW50c1swXS5sZW5ndGhcbiAgaWYoZCA8IDEpIHtcbiAgICByZXR1cm4gW11cbiAgfVxuXG4gIC8vU3BlY2lhbCBjYXNlOiAgRm9yIDFEIHdlIGNhbiBqdXN0IHNvcnQgdGhlIHBvaW50c1xuICBpZihkID09PSAxKSB7XG4gICAgcmV0dXJuIHRyaWFuZ3VsYXRlMUQobiwgcG9pbnRzLCBpbmNsdWRlUG9pbnRBdEluZmluaXR5KVxuICB9XG4gIFxuICAvL0xpZnQgcG9pbnRzLCBzb3J0XG4gIHZhciBsaWZ0ZWQgPSBuZXcgQXJyYXkobilcbiAgdmFyIHVwcGVyID0gMS4wXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIHZhciBwID0gcG9pbnRzW2ldXG4gICAgdmFyIHggPSBuZXcgQXJyYXkoZCsxKVxuICAgIHZhciBsID0gMC4wXG4gICAgZm9yKHZhciBqPTA7IGo8ZDsgKytqKSB7XG4gICAgICB2YXIgdiA9IHBbal1cbiAgICAgIHhbal0gPSB2XG4gICAgICBsICs9IHYgKiB2XG4gICAgfVxuICAgIHhbZF0gPSBsXG4gICAgbGlmdGVkW2ldID0gbmV3IExpZnRlZFBvaW50KHgsIGkpXG4gICAgdXBwZXIgPSBNYXRoLm1heChsLCB1cHBlcilcbiAgfVxuICB1bmlxKGxpZnRlZCwgY29tcGFyZUxpZnRlZClcbiAgXG4gIC8vRG91YmxlIHBvaW50c1xuICBuID0gbGlmdGVkLmxlbmd0aFxuXG4gIC8vQ3JlYXRlIG5ldyBsaXN0IG9mIHBvaW50c1xuICB2YXIgZHBvaW50cyA9IG5ldyBBcnJheShuICsgZCArIDEpXG4gIHZhciBkaW5kZXggPSBuZXcgQXJyYXkobiArIGQgKyAxKVxuXG4gIC8vQWRkIHN0ZWluZXIgcG9pbnRzIGF0IHRvcFxuICB2YXIgdSA9IChkKzEpICogKGQrMSkgKiB1cHBlclxuICB2YXIgeSA9IG5ldyBBcnJheShkKzEpXG4gIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICB5W2ldID0gMC4wXG4gIH1cbiAgeVtkXSA9IHVcblxuICBkcG9pbnRzWzBdID0geS5zbGljZSgpXG4gIGRpbmRleFswXSA9IC0xXG5cbiAgZm9yKHZhciBpPTA7IGk8PWQ7ICsraSkge1xuICAgIHZhciB4ID0geS5zbGljZSgpXG4gICAgeFtpXSA9IDFcbiAgICBkcG9pbnRzW2krMV0gPSB4XG4gICAgZGluZGV4W2krMV0gPSAtMVxuICB9XG5cbiAgLy9Db3B5IHJlc3Qgb2YgdGhlIHBvaW50cyBvdmVyXG4gIGZvcih2YXIgaT0wOyBpPG47ICsraSkge1xuICAgIHZhciBoID0gbGlmdGVkW2ldXG4gICAgZHBvaW50c1tpICsgZCArIDFdID0gaC5wb2ludFxuICAgIGRpbmRleFtpICsgZCArIDFdID0gIGguaW5kZXhcbiAgfVxuXG4gIC8vQ29uc3RydWN0IGNvbnZleCBodWxsXG4gIHZhciBodWxsID0gY2goZHBvaW50cywgZmFsc2UpXG4gIGlmKGluY2x1ZGVQb2ludEF0SW5maW5pdHkpIHtcbiAgICBodWxsID0gaHVsbC5maWx0ZXIoZnVuY3Rpb24oY2VsbCkge1xuICAgICAgdmFyIGNvdW50ID0gMFxuICAgICAgZm9yKHZhciBqPTA7IGo8PWQ7ICsraikge1xuICAgICAgICB2YXIgdiA9IGRpbmRleFtjZWxsW2pdXVxuICAgICAgICBpZih2IDwgMCkge1xuICAgICAgICAgIGlmKCsrY291bnQgPj0gMikge1xuICAgICAgICAgICAgcmV0dXJuIGZhbHNlXG4gICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIGNlbGxbal0gPSB2XG4gICAgICB9XG4gICAgICByZXR1cm4gdHJ1ZVxuICAgIH0pXG4gIH0gZWxzZSB7XG4gICAgaHVsbCA9IGh1bGwuZmlsdGVyKGZ1bmN0aW9uKGNlbGwpIHtcbiAgICAgIGZvcih2YXIgaT0wOyBpPD1kOyArK2kpIHtcbiAgICAgICAgdmFyIHYgPSBkaW5kZXhbY2VsbFtpXV1cbiAgICAgICAgaWYodiA8IDApIHtcbiAgICAgICAgICByZXR1cm4gZmFsc2VcbiAgICAgICAgfVxuICAgICAgICBjZWxsW2ldID0gdlxuICAgICAgfVxuICAgICAgcmV0dXJuIHRydWVcbiAgICB9KVxuICB9XG5cbiAgaWYoZCAmIDEpIHtcbiAgICBmb3IodmFyIGk9MDsgaTxodWxsLmxlbmd0aDsgKytpKSB7XG4gICAgICB2YXIgaCA9IGh1bGxbaV1cbiAgICAgIHZhciB4ID0gaFswXVxuICAgICAgaFswXSA9IGhbMV1cbiAgICAgIGhbMV0gPSB4XG4gICAgfVxuICB9XG5cbiAgcmV0dXJuIGh1bGxcbn0iLCJpZih0eXBlb2Ygd2luZG93LnBlcmZvcm1hbmNlID09PSBcIm9iamVjdFwiKSB7XG4gIGlmKHdpbmRvdy5wZXJmb3JtYW5jZS5ub3cpIHtcbiAgICBtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKCkgeyByZXR1cm4gd2luZG93LnBlcmZvcm1hbmNlLm5vdygpIH1cbiAgfSBlbHNlIGlmKHdpbmRvdy5wZXJmb3JtYW5jZS53ZWJraXROb3cpIHtcbiAgICBtb2R1bGUuZXhwb3J0cyA9IGZ1bmN0aW9uKCkgeyByZXR1cm4gd2luZG93LnBlcmZvcm1hbmNlLndlYmtpdE5vdygpIH1cbiAgfVxufSBlbHNlIGlmKERhdGUubm93KSB7XG4gIG1vZHVsZS5leHBvcnRzID0gRGF0ZS5ub3dcbn0gZWxzZSB7XG4gIG1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oKSB7IHJldHVybiAobmV3IERhdGUoKSkuZ2V0VGltZSgpIH1cbn1cbiIsIi8vQWRhcHRlZCBmcm9tIGhlcmU6IGh0dHBzOi8vZGV2ZWxvcGVyLm1vemlsbGEub3JnL2VuLVVTL2RvY3MvV2ViL1JlZmVyZW5jZS9FdmVudHMvd2hlZWw/cmVkaXJlY3Rsb2NhbGU9ZW4tVVMmcmVkaXJlY3RzbHVnPURPTSUyRk1vemlsbGFfZXZlbnRfcmVmZXJlbmNlJTJGd2hlZWxcblxudmFyIHByZWZpeCA9IFwiXCIsIF9hZGRFdmVudExpc3RlbmVyLCBvbndoZWVsLCBzdXBwb3J0O1xuXG4vLyBkZXRlY3QgZXZlbnQgbW9kZWxcbmlmICggd2luZG93LmFkZEV2ZW50TGlzdGVuZXIgKSB7XG4gIF9hZGRFdmVudExpc3RlbmVyID0gXCJhZGRFdmVudExpc3RlbmVyXCI7XG59IGVsc2Uge1xuICBfYWRkRXZlbnRMaXN0ZW5lciA9IFwiYXR0YWNoRXZlbnRcIjtcbiAgcHJlZml4ID0gXCJvblwiO1xufVxuXG4vLyBkZXRlY3QgYXZhaWxhYmxlIHdoZWVsIGV2ZW50XG5zdXBwb3J0ID0gXCJvbndoZWVsXCIgaW4gZG9jdW1lbnQuY3JlYXRlRWxlbWVudChcImRpdlwiKSA/IFwid2hlZWxcIiA6IC8vIE1vZGVybiBicm93c2VycyBzdXBwb3J0IFwid2hlZWxcIlxuICAgICAgICAgIGRvY3VtZW50Lm9ubW91c2V3aGVlbCAhPT0gdW5kZWZpbmVkID8gXCJtb3VzZXdoZWVsXCIgOiAvLyBXZWJraXQgYW5kIElFIHN1cHBvcnQgYXQgbGVhc3QgXCJtb3VzZXdoZWVsXCJcbiAgICAgICAgICBcIkRPTU1vdXNlU2Nyb2xsXCI7IC8vIGxldCdzIGFzc3VtZSB0aGF0IHJlbWFpbmluZyBicm93c2VycyBhcmUgb2xkZXIgRmlyZWZveFxuXG5mdW5jdGlvbiBfYWRkV2hlZWxMaXN0ZW5lciggZWxlbSwgZXZlbnROYW1lLCBjYWxsYmFjaywgdXNlQ2FwdHVyZSApIHtcbiAgZWxlbVsgX2FkZEV2ZW50TGlzdGVuZXIgXSggcHJlZml4ICsgZXZlbnROYW1lLCBzdXBwb3J0ID09IFwid2hlZWxcIiA/IGNhbGxiYWNrIDogZnVuY3Rpb24oIG9yaWdpbmFsRXZlbnQgKSB7XG4gICAgIW9yaWdpbmFsRXZlbnQgJiYgKCBvcmlnaW5hbEV2ZW50ID0gd2luZG93LmV2ZW50ICk7XG5cbiAgICAvLyBjcmVhdGUgYSBub3JtYWxpemVkIGV2ZW50IG9iamVjdFxuICAgIHZhciBldmVudCA9IHtcbiAgICAgIC8vIGtlZXAgYSByZWYgdG8gdGhlIG9yaWdpbmFsIGV2ZW50IG9iamVjdFxuICAgICAgb3JpZ2luYWxFdmVudDogb3JpZ2luYWxFdmVudCxcbiAgICAgIHRhcmdldDogb3JpZ2luYWxFdmVudC50YXJnZXQgfHwgb3JpZ2luYWxFdmVudC5zcmNFbGVtZW50LFxuICAgICAgdHlwZTogXCJ3aGVlbFwiLFxuICAgICAgZGVsdGFNb2RlOiBvcmlnaW5hbEV2ZW50LnR5cGUgPT0gXCJNb3pNb3VzZVBpeGVsU2Nyb2xsXCIgPyAwIDogMSxcbiAgICAgIGRlbHRhWDogMCxcbiAgICAgIGRlbGF0WjogMCxcbiAgICAgIHByZXZlbnREZWZhdWx0OiBmdW5jdGlvbigpIHtcbiAgICAgICAgb3JpZ2luYWxFdmVudC5wcmV2ZW50RGVmYXVsdCA/XG4gICAgICAgICAgb3JpZ2luYWxFdmVudC5wcmV2ZW50RGVmYXVsdCgpIDpcbiAgICAgICAgICBvcmlnaW5hbEV2ZW50LnJldHVyblZhbHVlID0gZmFsc2U7XG4gICAgICB9XG4gICAgfTtcbiAgICBcbiAgICAvLyBjYWxjdWxhdGUgZGVsdGFZIChhbmQgZGVsdGFYKSBhY2NvcmRpbmcgdG8gdGhlIGV2ZW50XG4gICAgaWYgKCBzdXBwb3J0ID09IFwibW91c2V3aGVlbFwiICkge1xuICAgICAgZXZlbnQuZGVsdGFZID0gLSAxLzQwICogb3JpZ2luYWxFdmVudC53aGVlbERlbHRhO1xuICAgICAgLy8gV2Via2l0IGFsc28gc3VwcG9ydCB3aGVlbERlbHRhWFxuICAgICAgb3JpZ2luYWxFdmVudC53aGVlbERlbHRhWCAmJiAoIGV2ZW50LmRlbHRhWCA9IC0gMS80MCAqIG9yaWdpbmFsRXZlbnQud2hlZWxEZWx0YVggKTtcbiAgICB9IGVsc2Uge1xuICAgICAgZXZlbnQuZGVsdGFZID0gb3JpZ2luYWxFdmVudC5kZXRhaWw7XG4gICAgfVxuXG4gICAgLy8gaXQncyB0aW1lIHRvIGZpcmUgdGhlIGNhbGxiYWNrXG4gICAgcmV0dXJuIGNhbGxiYWNrKCBldmVudCApO1xuICB9LCB1c2VDYXB0dXJlIHx8IGZhbHNlICk7XG59XG5cbm1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24oIGVsZW0sIGNhbGxiYWNrLCB1c2VDYXB0dXJlICkge1xuICBfYWRkV2hlZWxMaXN0ZW5lciggZWxlbSwgc3VwcG9ydCwgY2FsbGJhY2ssIHVzZUNhcHR1cmUgKTtcblxuICAvLyBoYW5kbGUgTW96TW91c2VQaXhlbFNjcm9sbCBpbiBvbGRlciBGaXJlZm94XG4gIGlmKCBzdXBwb3J0ID09IFwiRE9NTW91c2VTY3JvbGxcIiApIHtcbiAgICBfYWRkV2hlZWxMaXN0ZW5lciggZWxlbSwgXCJNb3pNb3VzZVBpeGVsU2Nyb2xsXCIsIGNhbGxiYWNrLCB1c2VDYXB0dXJlICk7XG4gIH1cbn07IiwiLy8gaHR0cDovL3BhdWxpcmlzaC5jb20vMjAxMS9yZXF1ZXN0YW5pbWF0aW9uZnJhbWUtZm9yLXNtYXJ0LWFuaW1hdGluZy9cbi8vIGh0dHA6Ly9teS5vcGVyYS5jb20vZW1vbGxlci9ibG9nLzIwMTEvMTIvMjAvcmVxdWVzdGFuaW1hdGlvbmZyYW1lLWZvci1zbWFydC1lci1hbmltYXRpbmdcbiBcbi8vIHJlcXVlc3RBbmltYXRpb25GcmFtZSBwb2x5ZmlsbCBieSBFcmlrIE3DtmxsZXIuIGZpeGVzIGZyb20gUGF1bCBJcmlzaCBhbmQgVGlubyBaaWpkZWxcbiBcbi8vIE1JVCBsaWNlbnNlXG52YXIgbGFzdFRpbWUgPSAwO1xudmFyIHZlbmRvcnMgPSBbJ21zJywgJ21veicsICd3ZWJraXQnLCAnbyddO1xuZm9yKHZhciB4ID0gMDsgeCA8IHZlbmRvcnMubGVuZ3RoICYmICF3aW5kb3cucmVxdWVzdEFuaW1hdGlvbkZyYW1lOyArK3gpIHtcbiAgICB3aW5kb3cucmVxdWVzdEFuaW1hdGlvbkZyYW1lID0gd2luZG93W3ZlbmRvcnNbeF0rJ1JlcXVlc3RBbmltYXRpb25GcmFtZSddO1xuICAgIHdpbmRvdy5jYW5jZWxBbmltYXRpb25GcmFtZSA9IHdpbmRvd1t2ZW5kb3JzW3hdKydDYW5jZWxBbmltYXRpb25GcmFtZSddIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHx8IHdpbmRvd1t2ZW5kb3JzW3hdKydDYW5jZWxSZXF1ZXN0QW5pbWF0aW9uRnJhbWUnXTtcbn1cblxuaWYgKCF3aW5kb3cucmVxdWVzdEFuaW1hdGlvbkZyYW1lKVxuICAgIHdpbmRvdy5yZXF1ZXN0QW5pbWF0aW9uRnJhbWUgPSBmdW5jdGlvbihjYWxsYmFjaywgZWxlbWVudCkge1xuICAgICAgICB2YXIgY3VyclRpbWUgPSBuZXcgRGF0ZSgpLmdldFRpbWUoKTtcbiAgICAgICAgdmFyIHRpbWVUb0NhbGwgPSBNYXRoLm1heCgwLCAxNiAtIChjdXJyVGltZSAtIGxhc3RUaW1lKSk7XG4gICAgICAgIHZhciBpZCA9IHdpbmRvdy5zZXRUaW1lb3V0KGZ1bmN0aW9uKCkgeyBjYWxsYmFjayhjdXJyVGltZSArIHRpbWVUb0NhbGwpOyB9LCBcbiAgICAgICAgICB0aW1lVG9DYWxsKTtcbiAgICAgICAgbGFzdFRpbWUgPSBjdXJyVGltZSArIHRpbWVUb0NhbGw7XG4gICAgICAgIHJldHVybiBpZDtcbiAgICB9O1xuXG5pZiAoIXdpbmRvdy5jYW5jZWxBbmltYXRpb25GcmFtZSlcbiAgICB3aW5kb3cuY2FuY2VsQW5pbWF0aW9uRnJhbWUgPSBmdW5jdGlvbihpZCkge1xuICAgICAgICBjbGVhclRpbWVvdXQoaWQpO1xuICAgIH07XG4iLCJcInVzZSBzdHJpY3RcIlxuXG5mdW5jdGlvbiBjb21waWxlU2VhcmNoKGZ1bmNOYW1lLCBwcmVkaWNhdGUsIHJldmVyc2VkLCBleHRyYUFyZ3MsIHVzZU5kYXJyYXksIGVhcmx5T3V0KSB7XG4gIHZhciBjb2RlID0gW1xuICAgIFwiZnVuY3Rpb24gXCIsIGZ1bmNOYW1lLCBcIihhLGwsaCxcIiwgZXh0cmFBcmdzLmpvaW4oXCIsXCIpLCAgXCIpe1wiLFxuZWFybHlPdXQgPyBcIlwiIDogXCJ2YXIgaT1cIiwgKHJldmVyc2VkID8gXCJsLTFcIiA6IFwiaCsxXCIpLFxuXCI7d2hpbGUobDw9aCl7XFxcbnZhciBtPShsK2gpPj4+MSx4PWFcIiwgdXNlTmRhcnJheSA/IFwiLmdldChtKVwiIDogXCJbbV1cIl1cbiAgaWYoZWFybHlPdXQpIHtcbiAgICBpZihwcmVkaWNhdGUuaW5kZXhPZihcImNcIikgPCAwKSB7XG4gICAgICBjb2RlLnB1c2goXCI7aWYoeD09PXkpe3JldHVybiBtfWVsc2UgaWYoeDw9eSl7XCIpXG4gICAgfSBlbHNlIHtcbiAgICAgIGNvZGUucHVzaChcIjt2YXIgcD1jKHgseSk7aWYocD09PTApe3JldHVybiBtfWVsc2UgaWYocDw9MCl7XCIpXG4gICAgfVxuICB9IGVsc2Uge1xuICAgIGNvZGUucHVzaChcIjtpZihcIiwgcHJlZGljYXRlLCBcIil7aT1tO1wiKVxuICB9XG4gIGlmKHJldmVyc2VkKSB7XG4gICAgY29kZS5wdXNoKFwibD1tKzF9ZWxzZXtoPW0tMX1cIilcbiAgfSBlbHNlIHtcbiAgICBjb2RlLnB1c2goXCJoPW0tMX1lbHNle2w9bSsxfVwiKVxuICB9XG4gIGNvZGUucHVzaChcIn1cIilcbiAgaWYoZWFybHlPdXQpIHtcbiAgICBjb2RlLnB1c2goXCJyZXR1cm4gLTF9O1wiKVxuICB9IGVsc2Uge1xuICAgIGNvZGUucHVzaChcInJldHVybiBpfTtcIilcbiAgfVxuICByZXR1cm4gY29kZS5qb2luKFwiXCIpXG59XG5cbmZ1bmN0aW9uIGNvbXBpbGVCb3VuZHNTZWFyY2gocHJlZGljYXRlLCByZXZlcnNlZCwgc3VmZml4LCBlYXJseU91dCkge1xuICB2YXIgcmVzdWx0ID0gbmV3IEZ1bmN0aW9uKFtcbiAgY29tcGlsZVNlYXJjaChcIkFcIiwgXCJ4XCIgKyBwcmVkaWNhdGUgKyBcInlcIiwgcmV2ZXJzZWQsIFtcInlcIl0sIGZhbHNlLCBlYXJseU91dCksXG4gIGNvbXBpbGVTZWFyY2goXCJCXCIsIFwieFwiICsgcHJlZGljYXRlICsgXCJ5XCIsIHJldmVyc2VkLCBbXCJ5XCJdLCB0cnVlLCBlYXJseU91dCksXG4gIGNvbXBpbGVTZWFyY2goXCJQXCIsIFwiYyh4LHkpXCIgKyBwcmVkaWNhdGUgKyBcIjBcIiwgcmV2ZXJzZWQsIFtcInlcIiwgXCJjXCJdLCBmYWxzZSwgZWFybHlPdXQpLFxuICBjb21waWxlU2VhcmNoKFwiUVwiLCBcImMoeCx5KVwiICsgcHJlZGljYXRlICsgXCIwXCIsIHJldmVyc2VkLCBbXCJ5XCIsIFwiY1wiXSwgdHJ1ZSwgZWFybHlPdXQpLFxuXCJmdW5jdGlvbiBkaXNwYXRjaEJzZWFyY2hcIiwgc3VmZml4LCBcIihhLHksYyxsLGgpe1xcXG5pZihhLnNoYXBlKXtcXFxuaWYodHlwZW9mKGMpPT09J2Z1bmN0aW9uJyl7XFxcbnJldHVybiBRKGEsKGw9PT11bmRlZmluZWQpPzA6bHwwLChoPT09dW5kZWZpbmVkKT9hLnNoYXBlWzBdLTE6aHwwLHksYylcXFxufWVsc2V7XFxcbnJldHVybiBCKGEsKGM9PT11bmRlZmluZWQpPzA6Y3wwLChsPT09dW5kZWZpbmVkKT9hLnNoYXBlWzBdLTE6bHwwLHkpXFxcbn19ZWxzZXtcXFxuaWYodHlwZW9mKGMpPT09J2Z1bmN0aW9uJyl7XFxcbnJldHVybiBQKGEsKGw9PT11bmRlZmluZWQpPzA6bHwwLChoPT09dW5kZWZpbmVkKT9hLmxlbmd0aC0xOmh8MCx5LGMpXFxcbn1lbHNle1xcXG5yZXR1cm4gQShhLChjPT09dW5kZWZpbmVkKT8wOmN8MCwobD09PXVuZGVmaW5lZCk/YS5sZW5ndGgtMTpsfDAseSlcXFxufX19XFxcbnJldHVybiBkaXNwYXRjaEJzZWFyY2hcIiwgc3VmZml4XS5qb2luKFwiXCIpKVxuICByZXR1cm4gcmVzdWx0KClcbn1cblxubW9kdWxlLmV4cG9ydHMgPSB7XG4gIGdlOiBjb21waWxlQm91bmRzU2VhcmNoKFwiPj1cIiwgZmFsc2UsIFwiR0VcIiksXG4gIGd0OiBjb21waWxlQm91bmRzU2VhcmNoKFwiPlwiLCBmYWxzZSwgXCJHVFwiKSxcbiAgbHQ6IGNvbXBpbGVCb3VuZHNTZWFyY2goXCI8XCIsIHRydWUsIFwiTFRcIiksXG4gIGxlOiBjb21waWxlQm91bmRzU2VhcmNoKFwiPD1cIiwgdHJ1ZSwgXCJMRVwiKSxcbiAgZXE6IGNvbXBpbGVCb3VuZHNTZWFyY2goXCItXCIsIHRydWUsIFwiRVFcIiwgdHJ1ZSlcbn1cbiIsIi8qIVxuICAqIGRvbXJlYWR5IChjKSBEdXN0aW4gRGlheiAyMDE0IC0gTGljZW5zZSBNSVRcbiAgKi9cbiFmdW5jdGlvbiAobmFtZSwgZGVmaW5pdGlvbikge1xuXG4gIGlmICh0eXBlb2YgbW9kdWxlICE9ICd1bmRlZmluZWQnKSBtb2R1bGUuZXhwb3J0cyA9IGRlZmluaXRpb24oKVxuICBlbHNlIGlmICh0eXBlb2YgZGVmaW5lID09ICdmdW5jdGlvbicgJiYgdHlwZW9mIGRlZmluZS5hbWQgPT0gJ29iamVjdCcpIGRlZmluZShkZWZpbml0aW9uKVxuICBlbHNlIHRoaXNbbmFtZV0gPSBkZWZpbml0aW9uKClcblxufSgnZG9tcmVhZHknLCBmdW5jdGlvbiAoKSB7XG5cbiAgdmFyIGZucyA9IFtdLCBsaXN0ZW5lclxuICAgICwgZG9jID0gZG9jdW1lbnRcbiAgICAsIGRvbUNvbnRlbnRMb2FkZWQgPSAnRE9NQ29udGVudExvYWRlZCdcbiAgICAsIGxvYWRlZCA9IC9ebG9hZGVkfF5jLy50ZXN0KGRvYy5yZWFkeVN0YXRlKVxuXG4gIGlmICghbG9hZGVkKVxuICBkb2MuYWRkRXZlbnRMaXN0ZW5lcihkb21Db250ZW50TG9hZGVkLCBsaXN0ZW5lciA9IGZ1bmN0aW9uICgpIHtcbiAgICBkb2MucmVtb3ZlRXZlbnRMaXN0ZW5lcihkb21Db250ZW50TG9hZGVkLCBsaXN0ZW5lcilcbiAgICBsb2FkZWQgPSAxXG4gICAgd2hpbGUgKGxpc3RlbmVyID0gZm5zLnNoaWZ0KCkpIGxpc3RlbmVyKClcbiAgfSlcblxuICByZXR1cm4gZnVuY3Rpb24gKGZuKSB7XG4gICAgbG9hZGVkID8gZm4oKSA6IGZucy5wdXNoKGZuKVxuICB9XG5cbn0pO1xuIiwiXCJ1c2Ugc3RyaWN0XCJcblxuZnVuY3Rpb24gaW52ZXJ0KGhhc2gpIHtcbiAgdmFyIHJlc3VsdCA9IHt9XG4gIGZvcih2YXIgaSBpbiBoYXNoKSB7XG4gICAgaWYoaGFzaC5oYXNPd25Qcm9wZXJ0eShpKSkge1xuICAgICAgcmVzdWx0W2hhc2hbaV1dID0gaVxuICAgIH1cbiAgfVxuICByZXR1cm4gcmVzdWx0XG59XG5cbm1vZHVsZS5leHBvcnRzID0gaW52ZXJ0IiwiXCJ1c2Ugc3RyaWN0XCJcblxuZnVuY3Rpb24gaW90YShuKSB7XG4gIHZhciByZXN1bHQgPSBuZXcgQXJyYXkobilcbiAgZm9yKHZhciBpPTA7IGk8bjsgKytpKSB7XG4gICAgcmVzdWx0W2ldID0gaVxuICB9XG4gIHJldHVybiByZXN1bHRcbn1cblxubW9kdWxlLmV4cG9ydHMgPSBpb3RhIiwiXCJ1c2Ugc3RyaWN0XCJcblxuZnVuY3Rpb24gdW5pcXVlX3ByZWQobGlzdCwgY29tcGFyZSkge1xuICB2YXIgcHRyID0gMVxuICAgICwgbGVuID0gbGlzdC5sZW5ndGhcbiAgICAsIGE9bGlzdFswXSwgYj1saXN0WzBdXG4gIGZvcih2YXIgaT0xOyBpPGxlbjsgKytpKSB7XG4gICAgYiA9IGFcbiAgICBhID0gbGlzdFtpXVxuICAgIGlmKGNvbXBhcmUoYSwgYikpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbGlzdFtwdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGxpc3QubGVuZ3RoID0gcHRyXG4gIHJldHVybiBsaXN0XG59XG5cbmZ1bmN0aW9uIHVuaXF1ZV9lcShsaXN0KSB7XG4gIHZhciBwdHIgPSAxXG4gICAgLCBsZW4gPSBsaXN0Lmxlbmd0aFxuICAgICwgYT1saXN0WzBdLCBiID0gbGlzdFswXVxuICBmb3IodmFyIGk9MTsgaTxsZW47ICsraSwgYj1hKSB7XG4gICAgYiA9IGFcbiAgICBhID0gbGlzdFtpXVxuICAgIGlmKGEgIT09IGIpIHtcbiAgICAgIGlmKGkgPT09IHB0cikge1xuICAgICAgICBwdHIrK1xuICAgICAgICBjb250aW51ZVxuICAgICAgfVxuICAgICAgbGlzdFtwdHIrK10gPSBhXG4gICAgfVxuICB9XG4gIGxpc3QubGVuZ3RoID0gcHRyXG4gIHJldHVybiBsaXN0XG59XG5cbmZ1bmN0aW9uIHVuaXF1ZShsaXN0LCBjb21wYXJlLCBzb3J0ZWQpIHtcbiAgaWYobGlzdC5sZW5ndGggPT09IDApIHtcbiAgICByZXR1cm4gbGlzdFxuICB9XG4gIGlmKGNvbXBhcmUpIHtcbiAgICBpZighc29ydGVkKSB7XG4gICAgICBsaXN0LnNvcnQoY29tcGFyZSlcbiAgICB9XG4gICAgcmV0dXJuIHVuaXF1ZV9wcmVkKGxpc3QsIGNvbXBhcmUpXG4gIH1cbiAgaWYoIXNvcnRlZCkge1xuICAgIGxpc3Quc29ydCgpXG4gIH1cbiAgcmV0dXJuIHVuaXF1ZV9lcShsaXN0KVxufVxuXG5tb2R1bGUuZXhwb3J0cyA9IHVuaXF1ZVxuIiwidmFyIHVhID0gdHlwZW9mIHdpbmRvdyAhPT0gJ3VuZGVmaW5lZCcgPyB3aW5kb3cubmF2aWdhdG9yLnVzZXJBZ2VudCA6ICcnXG4gICwgaXNPU1ggPSAvT1MgWC8udGVzdCh1YSlcbiAgLCBpc09wZXJhID0gL09wZXJhLy50ZXN0KHVhKVxuICAsIG1heWJlRmlyZWZveCA9ICEvbGlrZSBHZWNrby8udGVzdCh1YSkgJiYgIWlzT3BlcmFcblxudmFyIGksIG91dHB1dCA9IG1vZHVsZS5leHBvcnRzID0ge1xuICAwOiAgaXNPU1ggPyAnPG1lbnU+JyA6ICc8VU5LPidcbiwgMTogICc8bW91c2UgMT4nXG4sIDI6ICAnPG1vdXNlIDI+J1xuLCAzOiAgJzxicmVhaz4nXG4sIDQ6ICAnPG1vdXNlIDM+J1xuLCA1OiAgJzxtb3VzZSA0PidcbiwgNjogICc8bW91c2UgNT4nXG4sIDg6ICAnPGJhY2tzcGFjZT4nXG4sIDk6ICAnPHRhYj4nXG4sIDEyOiAnPGNsZWFyPidcbiwgMTM6ICc8ZW50ZXI+J1xuLCAxNjogJzxzaGlmdD4nXG4sIDE3OiAnPGNvbnRyb2w+J1xuLCAxODogJzxhbHQ+J1xuLCAxOTogJzxwYXVzZT4nXG4sIDIwOiAnPGNhcHMtbG9jaz4nXG4sIDIxOiAnPGltZS1oYW5ndWw+J1xuLCAyMzogJzxpbWUtanVuamE+J1xuLCAyNDogJzxpbWUtZmluYWw+J1xuLCAyNTogJzxpbWUta2Fuamk+J1xuLCAyNzogJzxlc2NhcGU+J1xuLCAyODogJzxpbWUtY29udmVydD4nXG4sIDI5OiAnPGltZS1ub25jb252ZXJ0PidcbiwgMzA6ICc8aW1lLWFjY2VwdD4nXG4sIDMxOiAnPGltZS1tb2RlLWNoYW5nZT4nXG4sIDI3OiAnPGVzY2FwZT4nXG4sIDMyOiAnPHNwYWNlPidcbiwgMzM6ICc8cGFnZS11cD4nXG4sIDM0OiAnPHBhZ2UtZG93bj4nXG4sIDM1OiAnPGVuZD4nXG4sIDM2OiAnPGhvbWU+J1xuLCAzNzogJzxsZWZ0PidcbiwgMzg6ICc8dXA+J1xuLCAzOTogJzxyaWdodD4nXG4sIDQwOiAnPGRvd24+J1xuLCA0MTogJzxzZWxlY3Q+J1xuLCA0MjogJzxwcmludD4nXG4sIDQzOiAnPGV4ZWN1dGU+J1xuLCA0NDogJzxzbmFwc2hvdD4nXG4sIDQ1OiAnPGluc2VydD4nXG4sIDQ2OiAnPGRlbGV0ZT4nXG4sIDQ3OiAnPGhlbHA+J1xuLCA5MTogJzxtZXRhPicgIC8vIG1ldGEtbGVmdCAtLSBubyBvbmUgaGFuZGxlcyBsZWZ0IGFuZCByaWdodCBwcm9wZXJseSwgc28gd2UgY29lcmNlIGludG8gb25lLlxuLCA5MjogJzxtZXRhPicgIC8vIG1ldGEtcmlnaHRcbiwgOTM6IGlzT1NYID8gJzxtZXRhPicgOiAnPG1lbnU+JyAgICAgIC8vIGNocm9tZSxvcGVyYSxzYWZhcmkgYWxsIHJlcG9ydCB0aGlzIGZvciBtZXRhLXJpZ2h0IChvc3ggbWJwKS5cbiwgOTU6ICc8c2xlZXA+J1xuLCAxMDY6ICc8bnVtLSo+J1xuLCAxMDc6ICc8bnVtLSs+J1xuLCAxMDg6ICc8bnVtLWVudGVyPidcbiwgMTA5OiAnPG51bS0tPidcbiwgMTEwOiAnPG51bS0uPidcbiwgMTExOiAnPG51bS0vPidcbiwgMTQ0OiAnPG51bS1sb2NrPidcbiwgMTQ1OiAnPHNjcm9sbC1sb2NrPidcbiwgMTYwOiAnPHNoaWZ0LWxlZnQ+J1xuLCAxNjE6ICc8c2hpZnQtcmlnaHQ+J1xuLCAxNjI6ICc8Y29udHJvbC1sZWZ0PidcbiwgMTYzOiAnPGNvbnRyb2wtcmlnaHQ+J1xuLCAxNjQ6ICc8YWx0LWxlZnQ+J1xuLCAxNjU6ICc8YWx0LXJpZ2h0PidcbiwgMTY2OiAnPGJyb3dzZXItYmFjaz4nXG4sIDE2NzogJzxicm93c2VyLWZvcndhcmQ+J1xuLCAxNjg6ICc8YnJvd3Nlci1yZWZyZXNoPidcbiwgMTY5OiAnPGJyb3dzZXItc3RvcD4nXG4sIDE3MDogJzxicm93c2VyLXNlYXJjaD4nXG4sIDE3MTogJzxicm93c2VyLWZhdm9yaXRlcz4nXG4sIDE3MjogJzxicm93c2VyLWhvbWU+J1xuXG4gIC8vIGZmL29zeCByZXBvcnRzICc8dm9sdW1lLW11dGU+JyBmb3IgJy0nXG4sIDE3MzogaXNPU1ggJiYgbWF5YmVGaXJlZm94ID8gJy0nIDogJzx2b2x1bWUtbXV0ZT4nXG4sIDE3NDogJzx2b2x1bWUtZG93bj4nXG4sIDE3NTogJzx2b2x1bWUtdXA+J1xuLCAxNzY6ICc8bmV4dC10cmFjaz4nXG4sIDE3NzogJzxwcmV2LXRyYWNrPidcbiwgMTc4OiAnPHN0b3A+J1xuLCAxNzk6ICc8cGxheS1wYXVzZT4nXG4sIDE4MDogJzxsYXVuY2gtbWFpbD4nXG4sIDE4MTogJzxsYXVuY2gtbWVkaWEtc2VsZWN0PidcbiwgMTgyOiAnPGxhdW5jaC1hcHAgMT4nXG4sIDE4MzogJzxsYXVuY2gtYXBwIDI+J1xuLCAxODY6ICc7J1xuLCAxODc6ICc9J1xuLCAxODg6ICcsJ1xuLCAxODk6ICctJ1xuLCAxOTA6ICcuJ1xuLCAxOTE6ICcvJ1xuLCAxOTI6ICdgJ1xuLCAyMTk6ICdbJ1xuLCAyMjA6ICdcXFxcJ1xuLCAyMjE6ICddJ1xuLCAyMjI6IFwiJ1wiXG4sIDIyMzogJzxtZXRhPidcbiwgMjI0OiAnPG1ldGE+JyAgICAgICAvLyBmaXJlZm94IHJlcG9ydHMgbWV0YSBoZXJlLlxuLCAyMjY6ICc8YWx0LWdyPidcbiwgMjI5OiAnPGltZS1wcm9jZXNzPidcbiwgMjMxOiBpc09wZXJhID8gJ2AnIDogJzx1bmljb2RlPidcbiwgMjQ2OiAnPGF0dGVudGlvbj4nXG4sIDI0NzogJzxjcnNlbD4nXG4sIDI0ODogJzxleHNlbD4nXG4sIDI0OTogJzxlcmFzZS1lb2Y+J1xuLCAyNTA6ICc8cGxheT4nXG4sIDI1MTogJzx6b29tPidcbiwgMjUyOiAnPG5vLW5hbWU+J1xuLCAyNTM6ICc8cGEtMT4nXG4sIDI1NDogJzxjbGVhcj4nXG59XG5cbmZvcihpID0gNTg7IGkgPCA2NTsgKytpKSB7XG4gIG91dHB1dFtpXSA9IFN0cmluZy5mcm9tQ2hhckNvZGUoaSlcbn1cblxuLy8gMC05XG5mb3IoaSA9IDQ4OyBpIDwgNTg7ICsraSkge1xuICBvdXRwdXRbaV0gPSAoaSAtIDQ4KSsnJ1xufVxuXG4vLyBBLVpcbmZvcihpID0gNjU7IGkgPCA5MTsgKytpKSB7XG4gIG91dHB1dFtpXSA9IFN0cmluZy5mcm9tQ2hhckNvZGUoaSlcbn1cblxuLy8gbnVtMC05XG5mb3IoaSA9IDk2OyBpIDwgMTA2OyArK2kpIHtcbiAgb3V0cHV0W2ldID0gJzxudW0tJysoaSAtIDk2KSsnPidcbn1cblxuLy8gRjEtRjI0XG5mb3IoaSA9IDExMjsgaSA8IDEzNjsgKytpKSB7XG4gIG91dHB1dFtpXSA9ICdGJysoaS0xMTEpXG59XG4iLCJcInVzZSBzdHJpY3RcIlxuXG52YXIgRXZlbnRFbWl0dGVyID0gcmVxdWlyZShcImV2ZW50c1wiKS5FdmVudEVtaXR0ZXJcbiAgLCB1dGlsICAgICAgICAgPSByZXF1aXJlKFwidXRpbFwiKVxuICAsIGRvbXJlYWR5ICAgICA9IHJlcXVpcmUoXCJkb21yZWFkeVwiKVxuICAsIHZrZXkgICAgICAgICA9IHJlcXVpcmUoXCJ2a2V5XCIpXG4gICwgaW52ZXJ0ICAgICAgID0gcmVxdWlyZShcImludmVydC1oYXNoXCIpXG4gICwgdW5pcSAgICAgICAgID0gcmVxdWlyZShcInVuaXFcIilcbiAgLCBic2VhcmNoICAgICAgPSByZXF1aXJlKFwiYmluYXJ5LXNlYXJjaC1ib3VuZHNcIilcbiAgLCBpb3RhICAgICAgICAgPSByZXF1aXJlKFwiaW90YS1hcnJheVwiKVxuICAsIG1pbiAgICAgICAgICA9IE1hdGgubWluXG5cbi8vQnJvd3NlciBjb21wYXRpYmlsaXR5IGhhY2tzXG5yZXF1aXJlKFwiLi9saWIvcmFmLXBvbHlmaWxsLmpzXCIpXG52YXIgYWRkTW91c2VXaGVlbCA9IHJlcXVpcmUoXCIuL2xpYi9tb3VzZXdoZWVsLXBvbHlmaWxsLmpzXCIpXG52YXIgaHJ0aW1lID0gcmVxdWlyZShcIi4vbGliL2hydGltZS1wb2x5ZmlsbC5qc1wiKVxuXG4vL1JlbW92ZSBhbmdsZSBicmFjZXMgYW5kIG90aGVyIHVzZWxlc3MgY3JhcFxudmFyIGZpbHRlcmVkX3ZrZXkgPSAoZnVuY3Rpb24oKSB7XG4gIHZhciByZXN1bHQgPSBuZXcgQXJyYXkoMjU2KVxuICAgICwgaSwgaiwga1xuICBmb3IoaT0wOyBpPDI1NjsgKytpKSB7XG4gICAgcmVzdWx0W2ldID0gXCJVTktcIlxuICB9XG4gIGZvcihpIGluIHZrZXkpIHtcbiAgICBrID0gdmtleVtpXVxuICAgIGlmKGsuY2hhckF0KDApID09PSAnPCcgJiYgay5jaGFyQXQoay5sZW5ndGgtMSkgPT09ICc+Jykge1xuICAgICAgayA9IGsuc3Vic3RyaW5nKDEsIGsubGVuZ3RoLTEpXG4gICAgfVxuICAgIGsgPSBrLnJlcGxhY2UoL1xccy9nLCBcIi1cIilcbiAgICByZXN1bHRbcGFyc2VJbnQoaSldID0ga1xuICB9XG4gIHJldHVybiByZXN1bHRcbn0pKClcblxuLy9Db21wdXRlIG1pbmltYWwgY29tbW9uIHNldCBvZiBrZXlib2FyZCBmdW5jdGlvbnNcbnZhciBrZXlOYW1lcyA9IHVuaXEoT2JqZWN0LmtleXMoaW52ZXJ0KGZpbHRlcmVkX3ZrZXkpKSlcblxuLy9UcmFuc2xhdGVzIGEgdmlydHVhbCBrZXljb2RlIHRvIGEgbm9ybWFsaXplZCBrZXljb2RlXG5mdW5jdGlvbiB2aXJ0dWFsS2V5Q29kZShrZXkpIHtcbiAgcmV0dXJuIGJzZWFyY2guZXEoa2V5TmFtZXMsIGtleSlcbn1cblxuLy9NYXBzIGEgcGh5c2ljYWwga2V5Y29kZSB0byBhIG5vcm1hbGl6ZWQga2V5Y29kZVxuZnVuY3Rpb24gcGh5c2ljYWxLZXlDb2RlKGtleSkge1xuICByZXR1cm4gdmlydHVhbEtleUNvZGUoZmlsdGVyZWRfdmtleVtrZXldKVxufVxuXG4vL0dhbWUgc2hlbGxcbmZ1bmN0aW9uIEdhbWVTaGVsbCgpIHtcbiAgRXZlbnRFbWl0dGVyLmNhbGwodGhpcylcbiAgdGhpcy5fY3VyS2V5U3RhdGUgID0gbmV3IEFycmF5KGtleU5hbWVzLmxlbmd0aClcbiAgdGhpcy5fcHJlc3NDb3VudCAgID0gbmV3IEFycmF5KGtleU5hbWVzLmxlbmd0aClcbiAgdGhpcy5fcmVsZWFzZUNvdW50ID0gbmV3IEFycmF5KGtleU5hbWVzLmxlbmd0aClcbiAgXG4gIHRoaXMuX3RpY2tJbnRlcnZhbCA9IG51bGxcbiAgdGhpcy5fcmFmSGFuZGxlID0gbnVsbFxuICB0aGlzLl90aWNrUmF0ZSA9IDBcbiAgdGhpcy5fbGFzdFRpY2sgPSBocnRpbWUoKVxuICB0aGlzLl9mcmFtZVRpbWUgPSAwLjBcbiAgdGhpcy5fcGF1c2VkID0gdHJ1ZVxuICB0aGlzLl93aWR0aCA9IDBcbiAgdGhpcy5faGVpZ2h0ID0gMFxuICBcbiAgdGhpcy5fd2FudEZ1bGxzY3JlZW4gPSBmYWxzZVxuICB0aGlzLl93YW50UG9pbnRlckxvY2sgPSBmYWxzZVxuICB0aGlzLl9mdWxsc2NyZWVuQWN0aXZlID0gZmFsc2VcbiAgdGhpcy5fcG9pbnRlckxvY2tBY3RpdmUgPSBmYWxzZVxuICBcbiAgdGhpcy5fcmVuZGVyID0gcmVuZGVyLmJpbmQodW5kZWZpbmVkLCB0aGlzKVxuICBcbiAgZm9yKHZhciBpPTA7IGk8a2V5TmFtZXMubGVuZ3RoOyArK2kpIHtcbiAgICB0aGlzLl9jdXJLZXlTdGF0ZVtpXSA9IGZhbHNlXG4gICAgdGhpcy5fcHJlc3NDb3VudFtpXSA9IHRoaXMuX3JlbGVhc2VDb3VudFtpXSA9IDBcbiAgfVxuICBcbiAgLy9QdWJsaWMgbWVtYmVyc1xuICB0aGlzLmVsZW1lbnQgPSBudWxsXG4gIHRoaXMuYmluZGluZ3MgPSB7fVxuICB0aGlzLmZyYW1lU2tpcCA9IDEwMC4wXG4gIHRoaXMudGlja0NvdW50ID0gMFxuICB0aGlzLmZyYW1lQ291bnQgPSAwXG4gIHRoaXMuc3RhcnRUaW1lID0gaHJ0aW1lKClcbiAgdGhpcy50aWNrVGltZSA9IHRoaXMuX3RpY2tSYXRlXG4gIHRoaXMuZnJhbWVUaW1lID0gMTAuMFxuICB0aGlzLnN0aWNreUZ1bGxzY3JlZW4gPSBmYWxzZVxuICB0aGlzLnN0aWNreVBvaW50TG9jayA9IGZhbHNlXG4gIFxuICAvL1Njcm9sbCBzdHVmZlxuICB0aGlzLnNjcm9sbCA9IFswLDAsMF1cbiAgICBcbiAgLy9Nb3VzZSBzdGF0ZVxuICB0aGlzLm1vdXNlWCA9IDBcbiAgdGhpcy5tb3VzZVkgPSAwXG4gIHRoaXMucHJldk1vdXNlWCA9IDBcbiAgdGhpcy5wcmV2TW91c2VZID0gMFxufVxuXG51dGlsLmluaGVyaXRzKEdhbWVTaGVsbCwgRXZlbnRFbWl0dGVyKVxuXG52YXIgcHJvdG8gPSBHYW1lU2hlbGwucHJvdG90eXBlXG5cbi8vQmluZCBrZXluYW1lc1xucHJvdG8ua2V5TmFtZXMgPSBrZXlOYW1lc1xuXG4vL0JpbmRzIGEgdmlydHVhbCBrZXlib2FyZCBldmVudCB0byBhIHBoeXNpY2FsIGtleVxucHJvdG8uYmluZCA9IGZ1bmN0aW9uKHZpcnR1YWxfa2V5KSB7XG4gIC8vTG9vayB1cCBwcmV2aW91cyBrZXkgYmluZGluZ3NcbiAgdmFyIGFyclxuICBpZih2aXJ0dWFsX2tleSBpbiB0aGlzLmJpbmRpbmdzKSB7XG4gICAgYXJyID0gdGhpcy5iaW5kaW5nc1t2aXJ0dWFsX2tleV1cbiAgfSBlbHNlIHtcbiAgICBhcnIgPSBbXVxuICB9XG4gIC8vQWRkIGtleXMgdG8gbGlzdFxuICB2YXIgcGh5c2ljYWxfa2V5XG4gIGZvcih2YXIgaT0xLCBuPWFyZ3VtZW50cy5sZW5ndGg7IGk8bjsgKytpKSB7XG4gICAgcGh5c2ljYWxfa2V5ID0gYXJndW1lbnRzW2ldXG4gICAgaWYodmlydHVhbEtleUNvZGUocGh5c2ljYWxfa2V5KSA+PSAwKSB7XG4gICAgICBhcnIucHVzaChwaHlzaWNhbF9rZXkpXG4gICAgfSBlbHNlIGlmKHBoeXNpY2FsX2tleSBpbiB0aGlzLmJpbmRpbmdzKSB7XG4gICAgICB2YXIga2V5YmluZHMgPSB0aGlzLmJpbmRpbmdzW3BoeXNpY2FsX2tleV1cbiAgICAgIGZvcih2YXIgaj0wOyBqPGtleWJpbmRzLmxlbmd0aDsgKytqKSB7XG4gICAgICAgIGFyci5wdXNoKGtleWJpbmRzW2pdKVxuICAgICAgfVxuICAgIH1cbiAgfVxuICAvL1JlbW92ZSBhbnkgZHVwbGljYXRlIGtleXNcbiAgYXJyID0gdW5pcShhcnIpXG4gIGlmKGFyci5sZW5ndGggPiAwKSB7XG4gICAgdGhpcy5iaW5kaW5nc1t2aXJ0dWFsX2tleV0gPSBhcnJcbiAgfVxuICB0aGlzLmVtaXQoJ2JpbmQnLCB2aXJ0dWFsX2tleSwgYXJyKVxufVxuXG4vL1VuYmluZHMgYSB2aXJ0dWFsIGtleWJvYXJkIGV2ZW50XG5wcm90by51bmJpbmQgPSBmdW5jdGlvbih2aXJ0dWFsX2tleSkge1xuICBpZih2aXJ0dWFsX2tleSBpbiB0aGlzLmJpbmRpbmdzKSB7XG4gICAgZGVsZXRlIHRoaXMuYmluZGluZ3NbdmlydHVhbF9rZXldXG4gIH1cbiAgdGhpcy5lbWl0KCd1bmJpbmQnLCB2aXJ0dWFsX2tleSlcbn1cblxuLy9DaGVja3MgaWYgYSBrZXkgaXMgc2V0IGluIGEgZ2l2ZW4gc3RhdGVcbmZ1bmN0aW9uIGxvb2t1cEtleShzdGF0ZSwgYmluZGluZ3MsIGtleSkge1xuICBpZihrZXkgaW4gYmluZGluZ3MpIHtcbiAgICB2YXIgYXJyID0gYmluZGluZ3Nba2V5XVxuICAgIGZvcih2YXIgaT0wLCBuPWFyci5sZW5ndGg7IGk8bjsgKytpKSB7XG4gICAgICBpZihzdGF0ZVt2aXJ0dWFsS2V5Q29kZShhcnJbaV0pXSkge1xuICAgICAgICByZXR1cm4gdHJ1ZVxuICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gZmFsc2VcbiAgfVxuICB2YXIga2MgPSB2aXJ0dWFsS2V5Q29kZShrZXkpXG4gIGlmKGtjID49IDApIHtcbiAgICByZXR1cm4gc3RhdGVba2NdXG4gIH1cbiAgcmV0dXJuIGZhbHNlXG59XG5cbi8vQ2hlY2tzIGlmIGEga2V5IGlzIHNldCBpbiBhIGdpdmVuIHN0YXRlXG5mdW5jdGlvbiBsb29rdXBDb3VudChzdGF0ZSwgYmluZGluZ3MsIGtleSkge1xuICBpZihrZXkgaW4gYmluZGluZ3MpIHtcbiAgICB2YXIgYXJyID0gYmluZGluZ3Nba2V5XSwgciA9IDBcbiAgICBmb3IodmFyIGk9MCwgbj1hcnIubGVuZ3RoOyBpPG47ICsraSkge1xuICAgICAgciArPSBzdGF0ZVt2aXJ0dWFsS2V5Q29kZShhcnJbaV0pXVxuICAgIH1cbiAgICByZXR1cm4gclxuICB9XG4gIHZhciBrYyA9IHZpcnR1YWxLZXlDb2RlKGtleSlcbiAgaWYoa2MgPj0gMCkge1xuICAgIHJldHVybiBzdGF0ZVtrY11cbiAgfVxuICByZXR1cm4gMFxufVxuXG4vL0NoZWNrcyBpZiBhIGtleSAoZWl0aGVyIHBoeXNpY2FsIG9yIHZpcnR1YWwpIGlzIGN1cnJlbnRseSBoZWxkIGRvd25cbnByb3RvLmRvd24gPSBmdW5jdGlvbihrZXkpIHtcbiAgcmV0dXJuIGxvb2t1cEtleSh0aGlzLl9jdXJLZXlTdGF0ZSwgdGhpcy5iaW5kaW5ncywga2V5KVxufVxuXG4vL0NoZWNrcyBpZiBhIGtleSB3YXMgZXZlciBkb3duXG5wcm90by53YXNEb3duID0gZnVuY3Rpb24oa2V5KSB7XG4gIHJldHVybiB0aGlzLmRvd24oa2V5KSB8fCAhIXRoaXMucHJlc3Moa2V5KVxufVxuXG4vL09wcG9zaXRlIG9mIGRvd25cbnByb3RvLnVwID0gZnVuY3Rpb24oa2V5KSB7XG4gIHJldHVybiAhdGhpcy5kb3duKGtleSlcbn1cblxuLy9DaGVja3MgaWYgYSBrZXkgd2FzIHJlbGVhc2VkIGR1cmluZyBwcmV2aW91cyBmcmFtZVxucHJvdG8ud2FzVXAgPSBmdW5jdGlvbihrZXkpIHtcbiAgcmV0dXJuIHRoaXMudXAoa2V5KSB8fCAhIXRoaXMucmVsZWFzZShrZXkpXG59XG5cbi8vUmV0dXJucyB0aGUgbnVtYmVyIG9mIHRpbWVzIGEga2V5IHdhcyBwcmVzc2VkIHNpbmNlIGxhc3QgdGlja1xucHJvdG8ucHJlc3MgPSBmdW5jdGlvbihrZXkpIHtcbiAgcmV0dXJuIGxvb2t1cENvdW50KHRoaXMuX3ByZXNzQ291bnQsIHRoaXMuYmluZGluZ3MsIGtleSlcbn1cblxuLy9SZXR1cm5zIHRoZSBudW1iZXIgb2YgdGltZXMgYSBrZXkgd2FzIHJlbGVhc2VkIHNpbmNlIGxhc3QgdGlja1xucHJvdG8ucmVsZWFzZSA9IGZ1bmN0aW9uKGtleSkge1xuICByZXR1cm4gbG9va3VwQ291bnQodGhpcy5fcmVsZWFzZUNvdW50LCB0aGlzLmJpbmRpbmdzLCBrZXkpXG59XG5cbi8vUGF1c2UvdW5wYXVzZSB0aGUgZ2FtZSBsb29wXG5PYmplY3QuZGVmaW5lUHJvcGVydHkocHJvdG8sIFwicGF1c2VkXCIsIHtcbiAgZ2V0OiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5fcGF1c2VkXG4gIH0sXG4gIHNldDogZnVuY3Rpb24oc3RhdGUpIHtcbiAgICB2YXIgbnMgPSAhIXN0YXRlXG4gICAgaWYobnMgIT09IHRoaXMuX3BhdXNlZCkge1xuICAgICAgaWYoIXRoaXMuX3BhdXNlZCkge1xuICAgICAgICB0aGlzLl9wYXVzZWQgPSB0cnVlXG4gICAgICAgIHRoaXMuX2ZyYW1lVGltZSA9IG1pbigxLjAsIChocnRpbWUoKSAtIHRoaXMuX2xhc3RUaWNrKSAvIHRoaXMuX3RpY2tSYXRlKVxuICAgICAgICBjbGVhckludGVydmFsKHRoaXMuX3RpY2tJbnRlcnZhbClcbiAgICAgICAgLy9jYW5jZWxBbmltYXRpb25GcmFtZSh0aGlzLl9yYWZIYW5kbGUpXG4gICAgICB9IGVsc2Uge1xuICAgICAgICB0aGlzLl9wYXVzZWQgPSBmYWxzZVxuICAgICAgICB0aGlzLl9sYXN0VGljayA9IGhydGltZSgpIC0gTWF0aC5mbG9vcih0aGlzLl9mcmFtZVRpbWUgKiB0aGlzLl90aWNrUmF0ZSlcbiAgICAgICAgdGhpcy5fdGlja0ludGVydmFsID0gc2V0SW50ZXJ2YWwodGljaywgdGhpcy5fdGlja1JhdGUsIHRoaXMpXG4gICAgICAgIHRoaXMuX3JhZkhhbmRsZSA9IHJlcXVlc3RBbmltYXRpb25GcmFtZSh0aGlzLl9yZW5kZXIpXG4gICAgICB9XG4gICAgfVxuICB9XG59KVxuXG4vL0Z1bGxzY3JlZW4gc3RhdGUgdG9nZ2xlXG5cbmZ1bmN0aW9uIHRyeUZ1bGxzY3JlZW4oc2hlbGwpIHtcbiAgLy9SZXF1ZXN0IGZ1bGwgc2NyZWVuXG4gIHZhciBlbGVtID0gc2hlbGwuZWxlbWVudFxuICBcbiAgaWYoc2hlbGwuX3dhbnRGdWxsc2NyZWVuICYmICFzaGVsbC5fZnVsbHNjcmVlbkFjdGl2ZSkge1xuICAgIHZhciBmcyA9IGVsZW0ucmVxdWVzdEZ1bGxzY3JlZW4gfHxcbiAgICAgICAgICAgICBlbGVtLnJlcXVlc3RGdWxsU2NyZWVuIHx8XG4gICAgICAgICAgICAgZWxlbS53ZWJraXRSZXF1ZXN0RnVsbHNjcmVlbiB8fFxuICAgICAgICAgICAgIGVsZW0ud2Via2l0UmVxdWVzdEZ1bGxTY3JlZW4gfHxcbiAgICAgICAgICAgICBlbGVtLm1velJlcXVlc3RGdWxsc2NyZWVuIHx8XG4gICAgICAgICAgICAgZWxlbS5tb3pSZXF1ZXN0RnVsbFNjcmVlbiB8fFxuICAgICAgICAgICAgIGZ1bmN0aW9uKCkge31cbiAgICBmcy5jYWxsKGVsZW0pXG4gIH1cbiAgaWYoc2hlbGwuX3dhbnRQb2ludGVyTG9jayAmJiAhc2hlbGwuX3BvaW50ZXJMb2NrQWN0aXZlKSB7XG4gICAgdmFyIHBsID0gIGVsZW0ucmVxdWVzdFBvaW50ZXJMb2NrIHx8XG4gICAgICAgICAgICAgIGVsZW0ud2Via2l0UmVxdWVzdFBvaW50ZXJMb2NrIHx8XG4gICAgICAgICAgICAgIGVsZW0ubW96UmVxdWVzdFBvaW50ZXJMb2NrIHx8XG4gICAgICAgICAgICAgIGVsZW0ubXNSZXF1ZXN0UG9pbnRlckxvY2sgfHxcbiAgICAgICAgICAgICAgZWxlbS5vUmVxdWVzdFBvaW50ZXJMb2NrIHx8XG4gICAgICAgICAgICAgIGZ1bmN0aW9uKCkge31cbiAgICBwbC5jYWxsKGVsZW0pXG4gIH1cbn1cblxudmFyIGNhbmNlbEZ1bGxzY3JlZW4gPSBkb2N1bWVudC5leGl0RnVsbHNjcmVlbiB8fFxuICAgICAgICAgICAgICAgICAgICAgICBkb2N1bWVudC5jYW5jZWxGdWxsc2NyZWVuIHx8ICAvL1doeSBjYW4gbm8gb25lIGFncmVlIG9uIHRoaXM/XG4gICAgICAgICAgICAgICAgICAgICAgIGRvY3VtZW50LmNhbmNlbEZ1bGxTY3JlZW4gfHxcbiAgICAgICAgICAgICAgICAgICAgICAgZG9jdW1lbnQud2Via2l0Q2FuY2VsRnVsbHNjcmVlbiB8fFxuICAgICAgICAgICAgICAgICAgICAgICBkb2N1bWVudC53ZWJraXRDYW5jZWxGdWxsU2NyZWVuIHx8XG4gICAgICAgICAgICAgICAgICAgICAgIGRvY3VtZW50Lm1vekNhbmNlbEZ1bGxzY3JlZW4gfHxcbiAgICAgICAgICAgICAgICAgICAgICAgZG9jdW1lbnQubW96Q2FuY2VsRnVsbFNjcmVlbiB8fFxuICAgICAgICAgICAgICAgICAgICAgICBmdW5jdGlvbigpe31cblxuT2JqZWN0LmRlZmluZVByb3BlcnR5KHByb3RvLCBcImZ1bGxzY3JlZW5cIiwge1xuICBnZXQ6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLl9mdWxsc2NyZWVuQWN0aXZlXG4gIH0sXG4gIHNldDogZnVuY3Rpb24oc3RhdGUpIHtcbiAgICB2YXIgbnMgPSAhIXN0YXRlXG4gICAgaWYoIW5zKSB7XG4gICAgICB0aGlzLl93YW50RnVsbHNjcmVlbiA9IGZhbHNlXG4gICAgICBjYW5jZWxGdWxsc2NyZWVuLmNhbGwoZG9jdW1lbnQpXG4gICAgfSBlbHNlIHtcbiAgICAgIHRoaXMuX3dhbnRGdWxsc2NyZWVuID0gdHJ1ZVxuICAgICAgdHJ5RnVsbHNjcmVlbih0aGlzKVxuICAgIH1cbiAgICByZXR1cm4gdGhpcy5fZnVsbHNjcmVlbkFjdGl2ZVxuICB9XG59KVxuXG5mdW5jdGlvbiBoYW5kbGVGdWxsc2NyZWVuKHNoZWxsKSB7XG4gIHNoZWxsLl9mdWxsc2NyZWVuQWN0aXZlID0gZG9jdW1lbnQuZnVsbHNjcmVlbiB8fFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGRvY3VtZW50Lm1vekZ1bGxTY3JlZW4gfHxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBkb2N1bWVudC53ZWJraXRJc0Z1bGxTY3JlZW4gfHxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBmYWxzZVxuICBpZighc2hlbGwuc3RpY2t5RnVsbHNjcmVlbiAmJiBzaGVsbC5fZnVsbHNjcmVlbkFjdGl2ZSkge1xuICAgIHNoZWxsLl93YW50RnVsbHNjcmVlbiA9IGZhbHNlXG4gIH1cbn1cblxuLy9Qb2ludGVyIGxvY2sgc3RhdGUgdG9nZ2xlXG52YXIgZXhpdFBvaW50ZXJMb2NrID0gZG9jdW1lbnQuZXhpdFBvaW50ZXJMb2NrIHx8XG4gICAgICAgICAgICAgICAgICAgICAgZG9jdW1lbnQud2Via2l0RXhpdFBvaW50ZXJMb2NrIHx8XG4gICAgICAgICAgICAgICAgICAgICAgZG9jdW1lbnQubW96RXhpdFBvaW50ZXJMb2NrIHx8XG4gICAgICAgICAgICAgICAgICAgICAgZnVuY3Rpb24oKSB7fVxuXG5PYmplY3QuZGVmaW5lUHJvcGVydHkocHJvdG8sIFwicG9pbnRlckxvY2tcIiwge1xuICBnZXQ6IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLl9wb2ludGVyTG9ja0FjdGl2ZVxuICB9LFxuICBzZXQ6IGZ1bmN0aW9uKHN0YXRlKSB7XG4gICAgdmFyIG5zID0gISFzdGF0ZVxuICAgIGlmKCFucykge1xuICAgICAgdGhpcy5fd2FudFBvaW50ZXJMb2NrID0gZmFsc2VcbiAgICAgIGV4aXRQb2ludGVyTG9jay5jYWxsKGRvY3VtZW50KVxuICAgIH0gZWxzZSB7XG4gICAgICB0aGlzLl93YW50UG9pbnRlckxvY2sgPSB0cnVlXG4gICAgICB0cnlGdWxsc2NyZWVuKHRoaXMpXG4gICAgfVxuICAgIHJldHVybiB0aGlzLl9wb2ludGVyTG9ja0FjdGl2ZVxuICB9XG59KVxuXG5mdW5jdGlvbiBoYW5kbGVQb2ludGVyTG9ja0NoYW5nZShzaGVsbCwgZXZlbnQpIHtcbiAgc2hlbGwuX3BvaW50ZXJMb2NrQWN0aXZlID0gc2hlbGwuZWxlbWVudCA9PT0gKFxuICAgICAgZG9jdW1lbnQucG9pbnRlckxvY2tFbGVtZW50IHx8XG4gICAgICBkb2N1bWVudC5tb3pQb2ludGVyTG9ja0VsZW1lbnQgfHxcbiAgICAgIGRvY3VtZW50LndlYmtpdFBvaW50ZXJMb2NrRWxlbWVudCB8fFxuICAgICAgbnVsbClcbiAgaWYoIXNoZWxsLnN0aWNreVBvaW50ZXJMb2NrICYmIHNoZWxsLl9wb2ludGVyTG9ja0FjdGl2ZSkge1xuICAgIHNoZWxsLl93YW50UG9pbnRlckxvY2sgPSBmYWxzZVxuICB9XG59XG5cbi8vV2lkdGggYW5kIGhlaWdodFxuT2JqZWN0LmRlZmluZVByb3BlcnR5KHByb3RvLCBcIndpZHRoXCIsIHtcbiAgZ2V0OiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5lbGVtZW50LmNsaWVudFdpZHRoXG4gIH1cbn0pXG5PYmplY3QuZGVmaW5lUHJvcGVydHkocHJvdG8sIFwiaGVpZ2h0XCIsIHtcbiAgZ2V0OiBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5lbGVtZW50LmNsaWVudEhlaWdodFxuICB9XG59KVxuXG4vL1NldCBrZXkgc3RhdGVcbmZ1bmN0aW9uIHNldEtleVN0YXRlKHNoZWxsLCBrZXksIHN0YXRlKSB7XG4gIHZhciBwcyA9IHNoZWxsLl9jdXJLZXlTdGF0ZVtrZXldXG4gIGlmKHBzICE9PSBzdGF0ZSkge1xuICAgIGlmKHN0YXRlKSB7XG4gICAgICBzaGVsbC5fcHJlc3NDb3VudFtrZXldKytcbiAgICB9IGVsc2Uge1xuICAgICAgc2hlbGwuX3JlbGVhc2VDb3VudFtrZXldKytcbiAgICB9XG4gICAgc2hlbGwuX2N1cktleVN0YXRlW2tleV0gPSBzdGF0ZVxuICB9XG59XG5cbi8vVGlja3MgdGhlIGdhbWUgc3RhdGUgb25lIHVwZGF0ZVxuZnVuY3Rpb24gdGljayhzaGVsbCkge1xuICB2YXIgc2tpcCA9IGhydGltZSgpICsgc2hlbGwuZnJhbWVTa2lwXG4gICAgLCBwQ291bnQgPSBzaGVsbC5fcHJlc3NDb3VudFxuICAgICwgckNvdW50ID0gc2hlbGwuX3JlbGVhc2VDb3VudFxuICAgICwgaSwgcywgdFxuICAgICwgdHIgPSBzaGVsbC5fdGlja1JhdGVcbiAgICAsIG4gPSBrZXlOYW1lcy5sZW5ndGhcbiAgd2hpbGUoIXNoZWxsLl9wYXVzZWQgJiZcbiAgICAgICAgaHJ0aW1lKCkgPj0gc2hlbGwuX2xhc3RUaWNrICsgdHIpIHtcbiAgICBcbiAgICAvL1NraXAgZnJhbWVzIGlmIHdlIGFyZSBvdmVyIGJ1ZGdldFxuICAgIGlmKGhydGltZSgpID4gc2tpcCkge1xuICAgICAgc2hlbGwuX2xhc3RUaWNrID0gaHJ0aW1lKCkgKyB0clxuICAgICAgcmV0dXJuXG4gICAgfVxuICAgIFxuICAgIC8vVGljayB0aGUgZ2FtZVxuICAgIHMgPSBocnRpbWUoKVxuICAgIHNoZWxsLmVtaXQoXCJ0aWNrXCIpXG4gICAgdCA9IGhydGltZSgpXG4gICAgc2hlbGwudGlja1RpbWUgPSB0IC0gc1xuICAgIFxuICAgIC8vVXBkYXRlIGNvdW50ZXJzIGFuZCB0aW1lXG4gICAgKytzaGVsbC50aWNrQ291bnRcbiAgICBzaGVsbC5fbGFzdFRpY2sgKz0gdHJcbiAgICBcbiAgICAvL1NoaWZ0IGlucHV0IHN0YXRlXG4gICAgZm9yKGk9MDsgaTxuOyArK2kpIHtcbiAgICAgIHBDb3VudFtpXSA9IHJDb3VudFtpXSA9IDBcbiAgICB9XG4gICAgaWYoc2hlbGwuX3BvaW50ZXJMb2NrQWN0aXZlKSB7XG4gICAgICBzaGVsbC5wcmV2TW91c2VYID0gc2hlbGwubW91c2VYID0gc2hlbGwud2lkdGg+PjFcbiAgICAgIHNoZWxsLnByZXZNb3VzZVkgPSBzaGVsbC5tb3VzZVkgPSBzaGVsbC5oZWlnaHQ+PjFcbiAgICB9IGVsc2Uge1xuICAgICAgc2hlbGwucHJldk1vdXNlWCA9IHNoZWxsLm1vdXNlWFxuICAgICAgc2hlbGwucHJldk1vdXNlWSA9IHNoZWxsLm1vdXNlWVxuICAgIH1cbiAgICBzaGVsbC5zY3JvbGxbMF0gPSBzaGVsbC5zY3JvbGxbMV0gPSBzaGVsbC5zY3JvbGxbMl0gPSAwXG4gIH1cbn1cblxuLy9SZW5kZXIgc3R1ZmZcbmZ1bmN0aW9uIHJlbmRlcihzaGVsbCkge1xuXG4gIC8vUmVxdWVzdCBuZXh0IGZyYW1lXG4gIHNoZWxsLl9yYWZIYW5kbGUgPSByZXF1ZXN0QW5pbWF0aW9uRnJhbWUoc2hlbGwuX3JlbmRlcilcblxuICAvL1RpY2sgdGhlIHNoZWxsXG4gIHRpY2soc2hlbGwpXG4gIFxuICAvL0NvbXB1dGUgZnJhbWUgdGltZVxuICB2YXIgZHRcbiAgaWYoc2hlbGwuX3BhdXNlZCkge1xuICAgIGR0ID0gc2hlbGwuX2ZyYW1lVGltZVxuICB9IGVsc2Uge1xuICAgIGR0ID0gbWluKDEuMCwgKGhydGltZSgpIC0gc2hlbGwuX2xhc3RUaWNrKSAvIHNoZWxsLl90aWNrUmF0ZSlcbiAgfVxuICBcbiAgLy9EcmF3IGEgZnJhbWVcbiAgKytzaGVsbC5mcmFtZUNvdW50XG4gIHZhciBzID0gaHJ0aW1lKClcbiAgc2hlbGwuZW1pdChcInJlbmRlclwiLCBkdClcbiAgdmFyIHQgPSBocnRpbWUoKVxuICBzaGVsbC5mcmFtZVRpbWUgPSB0IC0gc1xuICBcbn1cblxuZnVuY3Rpb24gaXNGb2N1c2VkKHNoZWxsKSB7XG4gIHJldHVybiAoZG9jdW1lbnQuYWN0aXZlRWxlbWVudCA9PT0gZG9jdW1lbnQuYm9keSkgfHxcbiAgICAgICAgIChkb2N1bWVudC5hY3RpdmVFbGVtZW50ID09PSBzaGVsbC5lbGVtZW50KVxufVxuXG4vL1NldCBrZXkgdXBcbmZ1bmN0aW9uIGhhbmRsZUtleVVwKHNoZWxsLCBldikge1xuICBldi5wcmV2ZW50RGVmYXVsdCgpXG4gIHZhciBrYyA9IHBoeXNpY2FsS2V5Q29kZShldi5rZXlDb2RlIHx8IGV2LmNoYXIgfHwgZXYud2hpY2ggfHwgZXYuY2hhckNvZGUpXG4gIGlmKGtjID49IDApIHtcbiAgICBzZXRLZXlTdGF0ZShzaGVsbCwga2MsIGZhbHNlKVxuICB9XG59XG5cbi8vU2V0IGtleSBkb3duXG5mdW5jdGlvbiBoYW5kbGVLZXlEb3duKHNoZWxsLCBldikge1xuICBpZighaXNGb2N1c2VkKHNoZWxsKSkge1xuICAgIHJldHVyblxuICB9XG4gIGlmKGV2Lm1ldGFLZXkpIHtcbiAgICAvL0hhY2s6IENsZWFyIGtleSBzdGF0ZSB3aGVuIG1ldGEgZ2V0cyBwcmVzc2VkIHRvIHByZXZlbnQga2V5cyBzdGlja2luZ1xuICAgIGhhbmRsZUJsdXIoc2hlbGwsIGV2KVxuICB9IGVsc2Uge1xuICAgIGV2LnByZXZlbnREZWZhdWx0KClcbiAgICB2YXIga2MgPSBwaHlzaWNhbEtleUNvZGUoZXYua2V5Q29kZSB8fCBldi5jaGFyIHx8IGV2LndoaWNoIHx8IGV2LmNoYXJDb2RlKVxuICAgIGlmKGtjID49IDApIHtcbiAgICAgIHNldEtleVN0YXRlKHNoZWxsLCBrYywgdHJ1ZSlcbiAgICB9XG4gIH1cbn1cblxuLy9Nb3VzZSBldmVudHMgYXJlIHJlYWxseSBhbm5veWluZ1xudmFyIG1vdXNlQ29kZXMgPSBpb3RhKDMyKS5tYXAoZnVuY3Rpb24obikge1xuICByZXR1cm4gdmlydHVhbEtleUNvZGUoXCJtb3VzZS1cIiArIChuKzEpKVxufSlcblxuZnVuY3Rpb24gc2V0TW91c2VCdXR0b25zKHNoZWxsLCBidXR0b25zKSB7XG4gIGZvcih2YXIgaT0wOyBpPDMyOyArK2kpIHtcbiAgICBzZXRLZXlTdGF0ZShzaGVsbCwgbW91c2VDb2Rlc1tpXSwgISEoYnV0dG9ucyAmICgxPDxpKSkpXG4gIH1cbn1cblxuZnVuY3Rpb24gaGFuZGxlTW91c2VNb3ZlKHNoZWxsLCBldikge1xuICBpZihzaGVsbC5fcG9pbnRlckxvY2tBY3RpdmUpIHtcbiAgICB2YXIgbW92ZW1lbnRYID0gZXYubW92ZW1lbnRYICAgICAgIHx8XG4gICAgICAgICAgICAgICAgICAgIGV2Lm1vek1vdmVtZW50WCAgICB8fFxuICAgICAgICAgICAgICAgICAgICBldi53ZWJraXRNb3ZlbWVudFggfHxcbiAgICAgICAgICAgICAgICAgICAgMCxcbiAgICAgICAgbW92ZW1lbnRZID0gZXYubW92ZW1lbnRZICAgICAgIHx8XG4gICAgICAgICAgICAgICAgICAgIGV2Lm1vek1vdmVtZW50WSAgICB8fFxuICAgICAgICAgICAgICAgICAgICBldi53ZWJraXRNb3ZlbWVudFkgfHxcbiAgICAgICAgICAgICAgICAgICAgMFxuICAgIHNoZWxsLm1vdXNlWCArPSBtb3ZlbWVudFhcbiAgICBzaGVsbC5tb3VzZVkgKz0gbW92ZW1lbnRZXG4gIH0gZWxzZSB7XG4gICAgc2hlbGwubW91c2VYID0gZXYuY2xpZW50WCAtIHNoZWxsLmVsZW1lbnQub2Zmc2V0TGVmdFxuICAgIHNoZWxsLm1vdXNlWSA9IGV2LmNsaWVudFkgLSBzaGVsbC5lbGVtZW50Lm9mZnNldFRvcFxuICB9XG4gIHJldHVybiBmYWxzZVxufVxuXG5mdW5jdGlvbiBoYW5kbGVNb3VzZURvd24oc2hlbGwsIGV2KSB7XG4gIHNldEtleVN0YXRlKHNoZWxsLCBtb3VzZUNvZGVzW2V2LmJ1dHRvbl0sIHRydWUpXG4gIHJldHVybiBmYWxzZVxufVxuXG5mdW5jdGlvbiBoYW5kbGVNb3VzZVVwKHNoZWxsLCBldikge1xuICBzZXRLZXlTdGF0ZShzaGVsbCwgbW91c2VDb2Rlc1tldi5idXR0b25dLCBmYWxzZSlcbiAgcmV0dXJuIGZhbHNlXG59XG5cbmZ1bmN0aW9uIGhhbmRsZU1vdXNlRW50ZXIoc2hlbGwsIGV2KSB7XG4gIGlmKHNoZWxsLl9wb2ludGVyTG9ja0FjdGl2ZSkge1xuICAgIHNoZWxsLnByZXZNb3VzZVggPSBzaGVsbC5tb3VzZVggPSBzaGVsbC53aWR0aD4+MVxuICAgIHNoZWxsLnByZXZNb3VzZVkgPSBzaGVsbC5tb3VzZVkgPSBzaGVsbC5oZWlnaHQ+PjFcbiAgfSBlbHNlIHtcbiAgICBzaGVsbC5wcmV2TW91c2VYID0gc2hlbGwubW91c2VYID0gZXYuY2xpZW50WCAtIHNoZWxsLmVsZW1lbnQub2Zmc2V0TGVmdFxuICAgIHNoZWxsLnByZXZNb3VzZVkgPSBzaGVsbC5tb3VzZVkgPSBldi5jbGllbnRZIC0gc2hlbGwuZWxlbWVudC5vZmZzZXRUb3BcbiAgfVxuICByZXR1cm4gZmFsc2Vcbn1cblxuZnVuY3Rpb24gaGFuZGxlTW91c2VMZWF2ZShzaGVsbCwgZXYpIHtcbiAgc2V0TW91c2VCdXR0b25zKHNoZWxsLCAwKVxuICByZXR1cm4gZmFsc2Vcbn1cblxuLy9IYW5kbGUgbW91c2Ugd2hlZWwgZXZlbnRzXG5mdW5jdGlvbiBoYW5kbGVNb3VzZVdoZWVsKHNoZWxsLCBldikge1xuICB2YXIgc2NhbGUgPSAxXG4gIHN3aXRjaChldi5kZWx0YU1vZGUpIHtcbiAgICBjYXNlIDA6IC8vUGl4ZWxcbiAgICAgIHNjYWxlID0gMVxuICAgIGJyZWFrXG4gICAgY2FzZSAxOiAvL0xpbmVcbiAgICAgIHNjYWxlID0gMTJcbiAgICBicmVha1xuICAgIGNhc2UgMjogLy9QYWdlXG4gICAgICAgc2NhbGUgPSBzaGVsbC5oZWlnaHRcbiAgICBicmVha1xuICB9XG4gIC8vQWRkIHNjcm9sbFxuICBzaGVsbC5zY3JvbGxbMF0gKz0gIGV2LmRlbHRhWCAqIHNjYWxlXG4gIHNoZWxsLnNjcm9sbFsxXSArPSAgZXYuZGVsdGFZICogc2NhbGVcbiAgc2hlbGwuc2Nyb2xsWzJdICs9IChldi5kZWx0YVogKiBzY2FsZSl8fDAuMFxuICByZXR1cm4gZmFsc2Vcbn1cblxuZnVuY3Rpb24gaGFuZGxlQ29udGV4TWVudShzaGVsbCwgZXYpIHtcbiAgcmV0dXJuIGZhbHNlXG59XG5cbmZ1bmN0aW9uIGhhbmRsZUJsdXIoc2hlbGwsIGV2KSB7XG4gIHZhciBuID0ga2V5TmFtZXMubGVuZ3RoXG4gICAgLCBjID0gc2hlbGwuX2N1cktleVN0YXRlXG4gICAgLCByID0gc2hlbGwuX3JlbGVhc2VDb3VudFxuICAgICwgaVxuICBmb3IoaT0wOyBpPG47ICsraSkge1xuICAgIGlmKGNbaV0pIHtcbiAgICAgICsrcltpXVxuICAgIH1cbiAgICBjW2ldID0gZmFsc2VcbiAgfVxuICByZXR1cm4gZmFsc2Vcbn1cblxuZnVuY3Rpb24gaGFuZGxlUmVzaXplRWxlbWVudChzaGVsbCwgZXYpIHtcbiAgdmFyIHcgPSBzaGVsbC5lbGVtZW50LmNsaWVudFdpZHRofDBcbiAgdmFyIGggPSBzaGVsbC5lbGVtZW50LmNsaWVudEhlaWdodHwwXG4gIGlmKCh3ICE9PSBzaGVsbC5fd2lkdGgpIHx8IChoICE9PSBzaGVsbC5faGVpZ2h0KSkge1xuICAgIHNoZWxsLl93aWR0aCA9IHdcbiAgICBzaGVsbC5faGVpZ2h0ID0gaFxuICAgIHNoZWxsLmVtaXQoXCJyZXNpemVcIiwgdywgaClcbiAgfVxufVxuXG5mdW5jdGlvbiBtYWtlRGVmYXVsdENvbnRhaW5lcigpIHtcbiAgdmFyIGNvbnRhaW5lciA9IGRvY3VtZW50LmNyZWF0ZUVsZW1lbnQoXCJkaXZcIilcbiAgY29udGFpbmVyLnRhYmluZGV4ID0gMVxuICBjb250YWluZXIuc3R5bGUucG9zaXRpb24gPSBcImFic29sdXRlXCJcbiAgY29udGFpbmVyLnN0eWxlLmxlZnQgPSBcIjBweFwiXG4gIGNvbnRhaW5lci5zdHlsZS5yaWdodCA9IFwiMHB4XCJcbiAgY29udGFpbmVyLnN0eWxlLnRvcCA9IFwiMHB4XCJcbiAgY29udGFpbmVyLnN0eWxlLmJvdHRvbSA9IFwiMHB4XCJcbiAgY29udGFpbmVyLnN0eWxlLmhlaWdodCA9IFwiMTAwJVwiXG4gIGNvbnRhaW5lci5zdHlsZS5vdmVyZmxvdyA9IFwiaGlkZGVuXCJcbiAgZG9jdW1lbnQuYm9keS5hcHBlbmRDaGlsZChjb250YWluZXIpXG4gIGRvY3VtZW50LmJvZHkuc3R5bGUub3ZlcmZsb3cgPSBcImhpZGRlblwiIC8vUHJldmVudCBib3VuY2VcbiAgZG9jdW1lbnQuYm9keS5zdHlsZS5oZWlnaHQgPSBcIjEwMCVcIlxuICByZXR1cm4gY29udGFpbmVyXG59XG5cbmZ1bmN0aW9uIGNyZWF0ZVNoZWxsKG9wdGlvbnMpIHtcbiAgb3B0aW9ucyA9IG9wdGlvbnMgfHwge31cbiAgXG4gIC8vQ2hlY2sgZnVsbHNjcmVlbiBhbmQgcG9pbnRlciBsb2NrIGZsYWdzXG4gIHZhciB1c2VGdWxsc2NyZWVuID0gISFvcHRpb25zLmZ1bGxzY3JlZW5cbiAgdmFyIHVzZVBvaW50ZXJMb2NrID0gdXNlRnVsbHNjcmVlblxuICBpZih0eXBlb2Ygb3B0aW9ucy5wb2ludGVyTG9jayAhPT0gdW5kZWZpbmVkKSB7XG4gICAgdXNlUG9pbnRlckxvY2sgPSAhIW9wdGlvbnMucG9pbnRlckxvY2tcbiAgfVxuICBcbiAgLy9DcmVhdGUgaW5pdGlhbCBzaGVsbFxuICB2YXIgc2hlbGwgPSBuZXcgR2FtZVNoZWxsKClcbiAgc2hlbGwuX3RpY2tSYXRlID0gb3B0aW9ucy50aWNrUmF0ZSB8fCAzMFxuICBzaGVsbC5mcmFtZVNraXAgPSBvcHRpb25zLmZyYW1lU2tpcCB8fCAoc2hlbGwuX3RpY2tSYXRlKzUpICogNVxuICBzaGVsbC5zdGlja3lGdWxsc2NyZWVuID0gISFvcHRpb25zLnN0aWNreUZ1bGxzY3JlZW4gfHwgISFvcHRpb25zLnN0aWNreVxuICBzaGVsbC5zdGlja3lQb2ludGVyTG9jayA9ICEhb3B0aW9ucy5zdGlja1BvaW50ZXJMb2NrIHx8ICFvcHRpb25zLnN0aWNreVxuICBcbiAgLy9TZXQgYmluZGluZ3NcbiAgaWYob3B0aW9ucy5iaW5kaW5ncykge1xuICAgIHNoZWxsLmJpbmRpbmdzID0gYmluZGluZ3NcbiAgfVxuICBcbiAgLy9XYWl0IGZvciBkb20gdG8gaW50aWFpbGl6ZVxuICBzZXRUaW1lb3V0KGZ1bmN0aW9uKCkgeyBkb21yZWFkeShmdW5jdGlvbiBpbml0R2FtZVNoZWxsKCkge1xuICAgIFxuICAgIC8vUmV0cmlldmUgZWxlbWVudFxuICAgIHZhciBlbGVtZW50ID0gb3B0aW9ucy5lbGVtZW50XG4gICAgaWYodHlwZW9mIGVsZW1lbnQgPT09IFwic3RyaW5nXCIpIHtcbiAgICAgIHZhciBlID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcihlbGVtZW50KVxuICAgICAgaWYoIWUpIHtcbiAgICAgICAgZSA9IGRvY3VtZW50LmdldEVsZW1lbnRCeUlkKGVsZW1lbnQpXG4gICAgICB9XG4gICAgICBpZighZSkge1xuICAgICAgICBlID0gZG9jdW1lbnQuZ2V0RWxlbWVudEJ5Q2xhc3MoZWxlbWVudClbMF1cbiAgICAgIH1cbiAgICAgIGlmKCFlKSB7XG4gICAgICAgIGUgPSBtYWtlRGVmYXVsdENvbnRhaW5lcigpXG4gICAgICB9XG4gICAgICBzaGVsbC5lbGVtZW50ID0gZVxuICAgIH0gZWxzZSBpZih0eXBlb2YgZWxlbWVudCA9PT0gXCJvYmplY3RcIiAmJiAhIWVsZW1lbnQpIHtcbiAgICAgIHNoZWxsLmVsZW1lbnQgPSBlbGVtZW50XG4gICAgfSBlbHNlIGlmKHR5cGVvZiBlbGVtZW50ID09PSBcImZ1bmN0aW9uXCIpIHtcbiAgICAgIHNoZWxsLmVsZW1lbnQgPSBlbGVtZW50KClcbiAgICB9IGVsc2Uge1xuICAgICAgc2hlbGwuZWxlbWVudCA9IG1ha2VEZWZhdWx0Q29udGFpbmVyKClcbiAgICB9XG4gICAgXG4gICAgLy9EaXNhYmxlIHVzZXItc2VsZWN0XG4gICAgaWYoc2hlbGwuZWxlbWVudC5zdHlsZSkge1xuICAgICAgc2hlbGwuZWxlbWVudC5zdHlsZVtcIi13ZWJraXQtdG91Y2gtY2FsbG91dFwiXSA9IFwibm9uZVwiXG4gICAgICBzaGVsbC5lbGVtZW50LnN0eWxlW1wiLXdlYmtpdC11c2VyLXNlbGVjdFwiXSA9IFwibm9uZVwiXG4gICAgICBzaGVsbC5lbGVtZW50LnN0eWxlW1wiLWtodG1sLXVzZXItc2VsZWN0XCJdID0gXCJub25lXCJcbiAgICAgIHNoZWxsLmVsZW1lbnQuc3R5bGVbXCItbW96LXVzZXItc2VsZWN0XCJdID0gXCJub25lXCJcbiAgICAgIHNoZWxsLmVsZW1lbnQuc3R5bGVbXCItbXMtdXNlci1zZWxlY3RcIl0gPSBcIm5vbmVcIlxuICAgICAgc2hlbGwuZWxlbWVudC5zdHlsZVtcInVzZXItc2VsZWN0XCJdID0gXCJub25lXCJcbiAgICB9XG4gICAgXG4gICAgLy9Ib29rIHJlc2l6ZSBoYW5kbGVyXG4gICAgc2hlbGwuX3dpZHRoID0gc2hlbGwuZWxlbWVudC5jbGllbnRXaWR0aFxuICAgIHNoZWxsLl9oZWlnaHQgPSBzaGVsbC5lbGVtZW50LmNsaWVudEhlaWdodFxuICAgIHZhciBoYW5kbGVSZXNpemUgPSBoYW5kbGVSZXNpemVFbGVtZW50LmJpbmQodW5kZWZpbmVkLCBzaGVsbClcbiAgICBpZih0eXBlb2YgTXV0YXRpb25PYnNlcnZlciAhPT0gXCJ1bmRlZmluZWRcIikge1xuICAgICAgdmFyIG9ic2VydmVyID0gbmV3IE11dGF0aW9uT2JzZXJ2ZXIoaGFuZGxlUmVzaXplKVxuICAgICAgb2JzZXJ2ZXIub2JzZXJ2ZShzaGVsbC5lbGVtZW50LCB7XG4gICAgICAgIGF0dHJpYnV0ZXM6IHRydWUsXG4gICAgICAgIHN1YnRyZWU6IHRydWVcbiAgICAgIH0pXG4gICAgfSBlbHNlIHtcbiAgICAgIHNoZWxsLmVsZW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcIkRPTVN1YnRyZWVNb2RpZmllZFwiLCBoYW5kbGVSZXNpemUsIGZhbHNlKVxuICAgIH1cbiAgICB3aW5kb3cuYWRkRXZlbnRMaXN0ZW5lcihcInJlc2l6ZVwiLCBoYW5kbGVSZXNpemUsIGZhbHNlKVxuICAgIFxuICAgIC8vSG9vayBrZXlib2FyZCBsaXN0ZW5lclxuICAgIHdpbmRvdy5hZGRFdmVudExpc3RlbmVyKFwia2V5ZG93blwiLCBoYW5kbGVLZXlEb3duLmJpbmQodW5kZWZpbmVkLCBzaGVsbCksIGZhbHNlKVxuICAgIHdpbmRvdy5hZGRFdmVudExpc3RlbmVyKFwia2V5dXBcIiwgaGFuZGxlS2V5VXAuYmluZCh1bmRlZmluZWQsIHNoZWxsKSwgZmFsc2UpXG4gICAgXG4gICAgLy9EaXNhYmxlIHJpZ2h0IGNsaWNrXG4gICAgc2hlbGwuZWxlbWVudC5vbmNvbnRleHRtZW51ID0gaGFuZGxlQ29udGV4TWVudS5iaW5kKHVuZGVmaW5lZCwgc2hlbGwpXG4gICAgXG4gICAgLy9Ib29rIG1vdXNlIGxpc3RlbmVyc1xuICAgIHNoZWxsLmVsZW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcIm1vdXNlZG93blwiLCBoYW5kbGVNb3VzZURvd24uYmluZCh1bmRlZmluZWQsIHNoZWxsKSwgZmFsc2UpXG4gICAgc2hlbGwuZWxlbWVudC5hZGRFdmVudExpc3RlbmVyKFwibW91c2V1cFwiLCBoYW5kbGVNb3VzZVVwLmJpbmQodW5kZWZpbmVkLCBzaGVsbCksIGZhbHNlKVxuICAgIHNoZWxsLmVsZW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcIm1vdXNlbW92ZVwiLCBoYW5kbGVNb3VzZU1vdmUuYmluZCh1bmRlZmluZWQsIHNoZWxsKSwgZmFsc2UpXG4gICAgc2hlbGwuZWxlbWVudC5hZGRFdmVudExpc3RlbmVyKFwibW91c2VlbnRlclwiLCBoYW5kbGVNb3VzZUVudGVyLmJpbmQodW5kZWZpbmVkLCBzaGVsbCksIGZhbHNlKVxuICAgIFxuICAgIC8vTW91c2UgbGVhdmVcbiAgICB2YXIgbGVhdmUgPSBoYW5kbGVNb3VzZUxlYXZlLmJpbmQodW5kZWZpbmVkLCBzaGVsbClcbiAgICBzaGVsbC5lbGVtZW50LmFkZEV2ZW50TGlzdGVuZXIoXCJtb3VzZWxlYXZlXCIsIGxlYXZlLCBmYWxzZSlcbiAgICBzaGVsbC5lbGVtZW50LmFkZEV2ZW50TGlzdGVuZXIoXCJtb3VzZW91dFwiLCBsZWF2ZSwgZmFsc2UpXG4gICAgd2luZG93LmFkZEV2ZW50TGlzdGVuZXIoXCJtb3VzZWxlYXZlXCIsIGxlYXZlLCBmYWxzZSlcbiAgICB3aW5kb3cuYWRkRXZlbnRMaXN0ZW5lcihcIm1vdXNlb3V0XCIsIGxlYXZlLCBmYWxzZSlcbiAgICBcbiAgICAvL0JsdXIgZXZlbnQgXG4gICAgdmFyIGJsdXIgPSBoYW5kbGVCbHVyLmJpbmQodW5kZWZpbmVkLCBzaGVsbClcbiAgICBzaGVsbC5lbGVtZW50LmFkZEV2ZW50TGlzdGVuZXIoXCJibHVyXCIsIGJsdXIsIGZhbHNlKVxuICAgIHNoZWxsLmVsZW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcImZvY3Vzb3V0XCIsIGJsdXIsIGZhbHNlKVxuICAgIHNoZWxsLmVsZW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcImZvY3VzXCIsIGJsdXIsIGZhbHNlKVxuICAgIHdpbmRvdy5hZGRFdmVudExpc3RlbmVyKFwiYmx1clwiLCBibHVyLCBmYWxzZSlcbiAgICB3aW5kb3cuYWRkRXZlbnRMaXN0ZW5lcihcImZvY3Vzb3V0XCIsIGJsdXIsIGZhbHNlKVxuICAgIHdpbmRvdy5hZGRFdmVudExpc3RlbmVyKFwiZm9jdXNcIiwgYmx1ciwgZmFsc2UpXG5cbiAgICAvL01vdXNlIHdoZWVsIGhhbmRsZXJcbiAgICBhZGRNb3VzZVdoZWVsKHNoZWxsLmVsZW1lbnQsIGhhbmRsZU1vdXNlV2hlZWwuYmluZCh1bmRlZmluZWQsIHNoZWxsKSwgZmFsc2UpXG5cbiAgICAvL0Z1bGxzY3JlZW4gaGFuZGxlclxuICAgIHZhciBmdWxsc2NyZWVuQ2hhbmdlID0gaGFuZGxlRnVsbHNjcmVlbi5iaW5kKHVuZGVmaW5lZCwgc2hlbGwpXG4gICAgZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcImZ1bGxzY3JlZW5jaGFuZ2VcIiwgZnVsbHNjcmVlbkNoYW5nZSwgZmFsc2UpXG4gICAgZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcIm1vemZ1bGxzY3JlZW5jaGFuZ2VcIiwgZnVsbHNjcmVlbkNoYW5nZSwgZmFsc2UpXG4gICAgZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcIndlYmtpdGZ1bGxzY3JlZW5jaGFuZ2VcIiwgZnVsbHNjcmVlbkNoYW5nZSwgZmFsc2UpXG5cbiAgICAvL1N0dXBpZCBmdWxsc2NyZWVuIGhhY2tcbiAgICBzaGVsbC5lbGVtZW50LmFkZEV2ZW50TGlzdGVuZXIoXCJjbGlja1wiLCB0cnlGdWxsc2NyZWVuLmJpbmQodW5kZWZpbmVkLCBzaGVsbCksIGZhbHNlKVxuXG4gICAgLy9Qb2ludGVyIGxvY2sgY2hhbmdlIGhhbmRsZXJcbiAgICB2YXIgcG9pbnRlckxvY2tDaGFuZ2UgPSBoYW5kbGVQb2ludGVyTG9ja0NoYW5nZS5iaW5kKHVuZGVmaW5lZCwgc2hlbGwpXG4gICAgZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcInBvaW50ZXJsb2NrY2hhbmdlXCIsIHBvaW50ZXJMb2NrQ2hhbmdlLCBmYWxzZSlcbiAgICBkb2N1bWVudC5hZGRFdmVudExpc3RlbmVyKFwibW96cG9pbnRlcmxvY2tjaGFuZ2VcIiwgcG9pbnRlckxvY2tDaGFuZ2UsIGZhbHNlKVxuICAgIGRvY3VtZW50LmFkZEV2ZW50TGlzdGVuZXIoXCJ3ZWJraXRwb2ludGVybG9ja2NoYW5nZVwiLCBwb2ludGVyTG9ja0NoYW5nZSwgZmFsc2UpXG4gICAgZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcInBvaW50ZXJsb2NrbG9zdFwiLCBwb2ludGVyTG9ja0NoYW5nZSwgZmFsc2UpXG4gICAgZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcIndlYmtpdHBvaW50ZXJsb2NrbG9zdFwiLCBwb2ludGVyTG9ja0NoYW5nZSwgZmFsc2UpXG4gICAgZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcihcIm1venBvaW50ZXJsb2NrbG9zdFwiLCBwb2ludGVyTG9ja0NoYW5nZSwgZmFsc2UpXG4gICAgXG4gICAgLy9VcGRhdGUgZmxhZ3NcbiAgICBzaGVsbC5mdWxsc2NyZWVuID0gdXNlRnVsbHNjcmVlblxuICAgIHNoZWxsLnBvaW50ZXJMb2NrID0gdXNlUG9pbnRlckxvY2tcbiAgXG4gICAgLy9EZWZhdWx0IG1vdXNlIGJ1dHRvbiBhbGlhc2VzXG4gICAgc2hlbGwuYmluZChcIm1vdXNlLWxlZnRcIiwgICBcIm1vdXNlLTFcIilcbiAgICBzaGVsbC5iaW5kKFwibW91c2UtcmlnaHRcIiwgIFwibW91c2UtM1wiKVxuICAgIHNoZWxsLmJpbmQoXCJtb3VzZS1taWRkbGVcIiwgXCJtb3VzZS0yXCIpXG4gICAgXG4gICAgLy9Jbml0aWFsaXplIHRpY2sgY291bnRlclxuICAgIHNoZWxsLl9sYXN0VGljayA9IGhydGltZSgpXG4gICAgc2hlbGwuc3RhcnRUaW1lID0gaHJ0aW1lKClcblxuICAgIC8vVW5wYXVzZSBzaGVsbFxuICAgIHNoZWxsLnBhdXNlZCA9IGZhbHNlXG4gICAgXG4gICAgLy9FbWl0IGluaXRpYWxpemUgZXZlbnRcbiAgICBzaGVsbC5lbWl0KFwiaW5pdFwiKVxuICB9KX0sIDApXG4gIFxuICByZXR1cm4gc2hlbGxcbn1cblxubW9kdWxlLmV4cG9ydHMgPSBjcmVhdGVTaGVsbFxuIiwiLyogKFRoZSBNSVQgTGljZW5zZSlcclxuICpcclxuICogQ29weXJpZ2h0IChjKSAyMDEyIEJyYW5kb24gQmVudmllIDxodHRwOi8vYmJlbnZpZS5jb20+XHJcbiAqXHJcbiAqIFBlcm1pc3Npb24gaXMgaGVyZWJ5IGdyYW50ZWQsIGZyZWUgb2YgY2hhcmdlLCB0byBhbnkgcGVyc29uIG9idGFpbmluZyBhIGNvcHkgb2YgdGhpcyBzb2Z0d2FyZSBhbmRcclxuICogYXNzb2NpYXRlZCBkb2N1bWVudGF0aW9uIGZpbGVzICh0aGUgJ1NvZnR3YXJlJyksIHRvIGRlYWwgaW4gdGhlIFNvZnR3YXJlIHdpdGhvdXQgcmVzdHJpY3Rpb24sXHJcbiAqIGluY2x1ZGluZyB3aXRob3V0IGxpbWl0YXRpb24gdGhlIHJpZ2h0cyB0byB1c2UsIGNvcHksIG1vZGlmeSwgbWVyZ2UsIHB1Ymxpc2gsIGRpc3RyaWJ1dGUsXHJcbiAqIHN1YmxpY2Vuc2UsIGFuZC9vciBzZWxsIGNvcGllcyBvZiB0aGUgU29mdHdhcmUsIGFuZCB0byBwZXJtaXQgcGVyc29ucyB0byB3aG9tIHRoZSBTb2Z0d2FyZSBpc1xyXG4gKiBmdXJuaXNoZWQgdG8gZG8gc28sIHN1YmplY3QgdG8gdGhlIGZvbGxvd2luZyBjb25kaXRpb25zOlxyXG4gKlxyXG4gKiBUaGUgYWJvdmUgY29weXJpZ2h0IG5vdGljZSBhbmQgdGhpcyBwZXJtaXNzaW9uIG5vdGljZSBzaGFsbCBiZSBpbmNsdWRlZCB3aXRoIGFsbCBjb3BpZXMgb3JcclxuICogc3Vic3RhbnRpYWwgcG9ydGlvbnMgb2YgdGhlIFNvZnR3YXJlLlxyXG4gKlxyXG4gKiBUSEUgU09GVFdBUkUgSVMgUFJPVklERUQgJ0FTIElTJywgV0lUSE9VVCBXQVJSQU5UWSBPRiBBTlkgS0lORCwgRVhQUkVTUyBPUiBJTVBMSUVELCBJTkNMVURJTkdcclxuICogQlVUIE5PVCBMSU1JVEVEIFRPIFRIRSBXQVJSQU5USUVTIE9GIE1FUkNIQU5UQUJJTElUWSwgRklUTkVTUyBGT1IgQSBQQVJUSUNVTEFSIFBVUlBPU0UgQU5EXHJcbiAqIE5PTklORlJJTkdFTUVOVC4gSU4gTk8gRVZFTlQgU0hBTEwgVEhFIEFVVEhPUlMgT1IgQ09QWVJJR0hUIEhPTERFUlMgQkUgTElBQkxFIEZPUiBBTlkgIENMQUlNLFxyXG4gKiBEQU1BR0VTIE9SIE9USEVSIExJQUJJTElUWSwgV0hFVEhFUiBJTiBBTiBBQ1RJT04gT0YgQ09OVFJBQ1QsIFRPUlQgT1IgT1RIRVJXSVNFLCBBUklTSU5HIEZST00sXHJcbiAqIE9VVCBPRiBPUiBJTiBDT05ORUNUSU9OIFdJVEggVEhFIFNPRlRXQVJFIE9SIFRIRSBVU0UgT1IgT1RIRVIgREVBTElOR1MgSU4gVEhFIFNPRlRXQVJFLlxyXG4gKi9cclxuXHJcbi8vIE9yaWdpbmFsIFdlYWtNYXAgaW1wbGVtZW50YXRpb24gYnkgR296YWxhIEAgaHR0cHM6Ly9naXN0LmdpdGh1Yi5jb20vMTI2OTk5MVxyXG4vLyBVcGRhdGVkIGFuZCBidWdmaXhlZCBieSBSYXlub3MgQCBodHRwczovL2dpc3QuZ2l0aHViLmNvbS8xNjM4MDU5XHJcbi8vIEV4cGFuZGVkIGJ5IEJlbnZpZSBAIGh0dHBzOi8vZ2l0aHViLmNvbS9CZW52aWUvaGFybW9ueS1jb2xsZWN0aW9uc1xyXG5cclxudm9pZCBmdW5jdGlvbihzdHJpbmdfLCBvYmplY3RfLCBmdW5jdGlvbl8sIHByb3RvdHlwZV8sIHRvU3RyaW5nXyxcclxuICAgICAgICAgICAgICBBcnJheSwgT2JqZWN0LCBGdW5jdGlvbiwgRlAsIGdsb2JhbCwgZXhwb3J0cywgdW5kZWZpbmVkXywgdW5kZWZpbmVkKXtcclxuXHJcbiAgdmFyIGdldFByb3BlcnRpZXMgPSBPYmplY3QuZ2V0T3duUHJvcGVydHlOYW1lcyxcclxuICAgICAgZXM1ID0gdHlwZW9mIGdldFByb3BlcnRpZXMgPT09IGZ1bmN0aW9uXyAmJiAhKHByb3RvdHlwZV8gaW4gZ2V0UHJvcGVydGllcyk7XHJcblxyXG4gIHZhciBjYWxsYmluZCA9IEZQLmJpbmRcclxuICAgID8gRlAuYmluZC5iaW5kKEZQLmNhbGwpXHJcbiAgICA6IChmdW5jdGlvbihjYWxsKXtcclxuICAgICAgICByZXR1cm4gZnVuY3Rpb24oZnVuYyl7XHJcbiAgICAgICAgICByZXR1cm4gZnVuY3Rpb24oKXtcclxuICAgICAgICAgICAgcmV0dXJuIGNhbGwuYXBwbHkoZnVuYywgYXJndW1lbnRzKTtcclxuICAgICAgICAgIH07XHJcbiAgICAgICAgfTtcclxuICAgICAgfShGUC5jYWxsKSk7XHJcblxyXG4gIHZhciBmdW5jdGlvblRvU3RyaW5nID0gY2FsbGJpbmQoRlBbdG9TdHJpbmdfXSksXHJcbiAgICAgIG9iamVjdFRvU3RyaW5nID0gY2FsbGJpbmQoe31bdG9TdHJpbmdfXSksXHJcbiAgICAgIG51bWJlclRvU3RyaW5nID0gY2FsbGJpbmQoLjAudG9TdHJpbmcpLFxyXG4gICAgICBjYWxsID0gY2FsbGJpbmQoRlAuY2FsbCksXHJcbiAgICAgIGFwcGx5ID0gY2FsbGJpbmQoRlAuYXBwbHkpLFxyXG4gICAgICBoYXNPd24gPSBjYWxsYmluZCh7fS5oYXNPd25Qcm9wZXJ0eSksXHJcbiAgICAgIHB1c2ggPSBjYWxsYmluZChbXS5wdXNoKSxcclxuICAgICAgc3BsaWNlID0gY2FsbGJpbmQoW10uc3BsaWNlKTtcclxuXHJcbiAgdmFyIG5hbWUgPSBmdW5jdGlvbihmdW5jKXtcclxuICAgIGlmICh0eXBlb2YgZnVuYyAhPT0gZnVuY3Rpb25fKVxyXG4gICAgICByZXR1cm4gJyc7XHJcbiAgICBlbHNlIGlmICgnbmFtZScgaW4gZnVuYylcclxuICAgICAgcmV0dXJuIGZ1bmMubmFtZTtcclxuXHJcbiAgICByZXR1cm4gZnVuY3Rpb25Ub1N0cmluZyhmdW5jKS5tYXRjaCgvXlxcbj9mdW5jdGlvblxccz8oXFx3Kik/Xz9cXCgvKVsxXTtcclxuICB9O1xyXG5cclxuICB2YXIgY3JlYXRlID0gZXM1XHJcbiAgICA/IE9iamVjdC5jcmVhdGVcclxuICAgIDogZnVuY3Rpb24ocHJvdG8sIGRlc2NzKXtcclxuICAgICAgICB2YXIgQ3RvciA9IGZ1bmN0aW9uKCl7fTtcclxuICAgICAgICBDdG9yW3Byb3RvdHlwZV9dID0gT2JqZWN0KHByb3RvKTtcclxuICAgICAgICB2YXIgb2JqZWN0ID0gbmV3IEN0b3I7XHJcblxyXG4gICAgICAgIGlmIChkZXNjcylcclxuICAgICAgICAgIGZvciAodmFyIGtleSBpbiBkZXNjcylcclxuICAgICAgICAgICAgZGVmaW5lUHJvcGVydHkob2JqZWN0LCBrZXksIGRlc2NzW2tdKTtcclxuXHJcbiAgICAgICAgcmV0dXJuIG9iamVjdDtcclxuICAgICAgfTtcclxuXHJcblxyXG4gIGZ1bmN0aW9uIEhhc2goKXt9XHJcblxyXG4gIGlmIChlczUpIHtcclxuICAgIHZvaWQgZnVuY3Rpb24oT2JqZWN0Q3JlYXRlKXtcclxuICAgICAgSGFzaC5wcm90b3R5cGUgPSBPYmplY3RDcmVhdGUobnVsbCk7XHJcbiAgICAgIGZ1bmN0aW9uIGluaGVyaXQob2JqKXtcclxuICAgICAgICByZXR1cm4gT2JqZWN0Q3JlYXRlKG9iaik7XHJcbiAgICAgIH1cclxuICAgICAgSGFzaC5pbmhlcml0ID0gaW5oZXJpdDtcclxuICAgIH0oT2JqZWN0LmNyZWF0ZSk7XHJcbiAgfSBlbHNlIHtcclxuICAgIHZvaWQgZnVuY3Rpb24oRil7XHJcbiAgICAgIHZhciBpZnJhbWUgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50KCdpZnJhbWUnKTtcclxuICAgICAgaWZyYW1lLnN0eWxlLmRpc3BsYXkgPSAnbm9uZSc7XHJcbiAgICAgIGRvY3VtZW50LmJvZHkuYXBwZW5kQ2hpbGQoaWZyYW1lKTtcclxuICAgICAgaWZyYW1lLnNyYyA9ICdqYXZhc2NyaXB0OidcclxuICAgICAgSGFzaC5wcm90b3R5cGUgPSBpZnJhbWUuY29udGVudFdpbmRvdy5PYmplY3QucHJvdG90eXBlO1xyXG4gICAgICBkb2N1bWVudC5ib2R5LnJlbW92ZUNoaWxkKGlmcmFtZSk7XHJcbiAgICAgIGlmcmFtZSA9IG51bGw7XHJcblxyXG4gICAgICB2YXIgcHJvcHMgPSBbJ2NvbnN0cnVjdG9yJywgJ2hhc093blByb3BlcnR5JywgJ3Byb3BlcnR5SXNFbnVtZXJhYmxlJyxcclxuICAgICAgICAgICAgICAgICAgICdpc1Byb3RveXBlT2YnLCAndG9Mb2NhbGVTdHJpbmcnLCAndG9TdHJpbmcnLCAndmFsdWVPZiddO1xyXG5cclxuICAgICAgZm9yICh2YXIgaT0wOyBpIDwgcHJvcHMubGVuZ3RoOyBpKyspXHJcbiAgICAgICAgZGVsZXRlIEhhc2gucHJvdG90eXBlW3Byb3BzW2ldXTtcclxuXHJcbiAgICAgIGZ1bmN0aW9uIGluaGVyaXQob2JqKXtcclxuICAgICAgICBGLnByb3RvdHlwZSA9IG9iajtcclxuICAgICAgICBvYmogPSBuZXcgRjtcclxuICAgICAgICBGLnByb3RvdHlwZSA9IG51bGw7XHJcbiAgICAgICAgcmV0dXJuIG9iajtcclxuICAgICAgfVxyXG5cclxuICAgICAgSGFzaC5pbmhlcml0ID0gaW5oZXJpdDtcclxuICAgIH0oZnVuY3Rpb24oKXt9KTtcclxuICB9XHJcblxyXG4gIHZhciBkZWZpbmVQcm9wZXJ0eSA9IGVzNVxyXG4gICAgPyBPYmplY3QuZGVmaW5lUHJvcGVydHlcclxuICAgIDogZnVuY3Rpb24ob2JqZWN0LCBrZXksIGRlc2MpIHtcclxuICAgICAgICBvYmplY3Rba2V5XSA9IGRlc2MudmFsdWU7XHJcbiAgICAgICAgcmV0dXJuIG9iamVjdDtcclxuICAgICAgfTtcclxuXHJcbiAgdmFyIGRlZmluZSA9IGZ1bmN0aW9uKG9iamVjdCwga2V5LCB2YWx1ZSl7XHJcbiAgICBpZiAodHlwZW9mIGtleSA9PT0gZnVuY3Rpb25fKSB7XHJcbiAgICAgIHZhbHVlID0ga2V5O1xyXG4gICAgICBrZXkgPSBuYW1lKHZhbHVlKS5yZXBsYWNlKC9fJC8sICcnKTtcclxuICAgIH1cclxuXHJcbiAgICByZXR1cm4gZGVmaW5lUHJvcGVydHkob2JqZWN0LCBrZXksIHsgY29uZmlndXJhYmxlOiB0cnVlLCB3cml0YWJsZTogdHJ1ZSwgdmFsdWU6IHZhbHVlIH0pO1xyXG4gIH07XHJcblxyXG4gIHZhciBpc0FycmF5ID0gZXM1XHJcbiAgICA/IChmdW5jdGlvbihpc0FycmF5KXtcclxuICAgICAgICByZXR1cm4gZnVuY3Rpb24obyl7XHJcbiAgICAgICAgICByZXR1cm4gaXNBcnJheShvKSB8fCBvIGluc3RhbmNlb2YgQXJyYXk7XHJcbiAgICAgICAgfTtcclxuICAgICAgfSkoQXJyYXkuaXNBcnJheSlcclxuICAgIDogZnVuY3Rpb24obyl7XHJcbiAgICAgICAgcmV0dXJuIG8gaW5zdGFuY2VvZiBBcnJheSB8fCBvYmplY3RUb1N0cmluZyhvKSA9PT0gJ1tvYmplY3QgQXJyYXldJztcclxuICAgICAgfTtcclxuXHJcbiAgLy8gIyMjIyMjIyMjIyMjXHJcbiAgLy8gIyMjIERhdGEgIyMjXHJcbiAgLy8gIyMjIyMjIyMjIyMjXHJcblxyXG4gIHZhciBidWlsdGluV2Vha01hcCA9ICdXZWFrTWFwJyBpbiBnbG9iYWw7XHJcblxyXG4gIHZhciBNYXBEYXRhID0gYnVpbHRpbldlYWtNYXBcclxuICAgID8gKGZ1bmN0aW9uKCl7XHJcbiAgICAgIHZhciBCdWlsdGluV2Vha01hcCA9IGdsb2JhbC5XZWFrTWFwLFxyXG4gICAgICAgICAgd21nZXQgPSBjYWxsYmluZChCdWlsdGluV2Vha01hcFtwcm90b3R5cGVfXS5nZXQpLFxyXG4gICAgICAgICAgd21zZXQgPSBjYWxsYmluZChCdWlsdGluV2Vha01hcFtwcm90b3R5cGVfXS5zZXQpLFxyXG4gICAgICAgICAgd21oYXMgPSBjYWxsYmluZChCdWlsdGluV2Vha01hcFtwcm90b3R5cGVfXS5oYXMpO1xyXG5cclxuICAgICAgZnVuY3Rpb24gTWFwRGF0YShuYW1lKXtcclxuICAgICAgICB2YXIgbWFwID0gbmV3IEJ1aWx0aW5XZWFrTWFwO1xyXG5cclxuICAgICAgICB0aGlzLmdldCA9IGZ1bmN0aW9uKG8pe1xyXG4gICAgICAgICAgcmV0dXJuIHdtZ2V0KG1hcCwgbyk7XHJcbiAgICAgICAgfTtcclxuICAgICAgICB0aGlzLnNldCA9IGZ1bmN0aW9uKG8sIHYpe1xyXG4gICAgICAgICAgd21zZXQobWFwLCBvLCB2KTtcclxuICAgICAgICB9O1xyXG5cclxuICAgICAgICBpZiAobmFtZSkge1xyXG4gICAgICAgICAgdGhpcy53cmFwID0gZnVuY3Rpb24obywgdil7XHJcbiAgICAgICAgICAgIGlmICh3bWhhcyhtYXAsIG8pKVxyXG4gICAgICAgICAgICAgIHRocm93IG5ldyBUeXBlRXJyb3IoXCJPYmplY3QgaXMgYWxyZWFkeSBhIFwiICsgbmFtZSk7XHJcbiAgICAgICAgICAgIHdtc2V0KG1hcCwgbywgdik7XHJcbiAgICAgICAgICB9O1xyXG4gICAgICAgICAgdGhpcy51bndyYXAgPSBmdW5jdGlvbihvKXtcclxuICAgICAgICAgICAgdmFyIHN0b3JhZ2UgPSB3bWdldChtYXAsIG8pO1xyXG4gICAgICAgICAgICBpZiAoIXN0b3JhZ2UpXHJcbiAgICAgICAgICAgICAgdGhyb3cgbmV3IFR5cGVFcnJvcihuYW1lICsgXCIgaXMgbm90IGdlbmVyaWNcIik7XHJcbiAgICAgICAgICAgIHJldHVybiBzdG9yYWdlO1xyXG4gICAgICAgICAgfTtcclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIHJldHVybiBNYXBEYXRhO1xyXG4gICAgfSkoKVxyXG4gICAgOiAoZnVuY3Rpb24oKXtcclxuICAgICAgdmFyIGxvY2tlciA9ICdyZXR1cm4gZnVuY3Rpb24oayl7aWYoaz09PXMpcmV0dXJuIGx9JyxcclxuICAgICAgICAgIHJhbmRvbSA9IE1hdGgucmFuZG9tLFxyXG4gICAgICAgICAgdWlkcyA9IG5ldyBIYXNoLFxyXG4gICAgICAgICAgc2xpY2UgPSBjYWxsYmluZCgnJy5zbGljZSksXHJcbiAgICAgICAgICBpbmRleE9mID0gY2FsbGJpbmQoW10uaW5kZXhPZik7XHJcblxyXG4gICAgICB2YXIgY3JlYXRlVUlEID0gZnVuY3Rpb24oKXtcclxuICAgICAgICB2YXIga2V5ID0gc2xpY2UobnVtYmVyVG9TdHJpbmcocmFuZG9tKCksIDM2KSwgMik7XHJcbiAgICAgICAgcmV0dXJuIGtleSBpbiB1aWRzID8gY3JlYXRlVUlEKCkgOiB1aWRzW2tleV0gPSBrZXk7XHJcbiAgICAgIH07XHJcblxyXG4gICAgICB2YXIgZ2xvYmFsSUQgPSBjcmVhdGVVSUQoKTtcclxuXHJcbiAgICAgIC8vIGNvbW1vbiBwZXItb2JqZWN0IHN0b3JhZ2UgYXJlYSBtYWRlIHZpc2libGUgYnkgcGF0Y2hpbmcgZ2V0T3duUHJvcGVydHlOYW1lcydcclxuICAgICAgZnVuY3Rpb24gZ2V0T3duUHJvcGVydHlOYW1lcyhvYmope1xyXG4gICAgICAgIHZhciBwcm9wcyA9IGdldFByb3BlcnRpZXMob2JqKTtcclxuICAgICAgICBpZiAoaGFzT3duKG9iaiwgZ2xvYmFsSUQpKVxyXG4gICAgICAgICAgc3BsaWNlKHByb3BzLCBpbmRleE9mKHByb3BzLCBnbG9iYWxJRCksIDEpO1xyXG4gICAgICAgIHJldHVybiBwcm9wcztcclxuICAgICAgfVxyXG5cclxuICAgICAgaWYgKGVzNSkge1xyXG4gICAgICAgIC8vIGNoZWNrIGZvciB0aGUgcmFuZG9tIGtleSBvbiBhbiBvYmplY3QsIGNyZWF0ZSBuZXcgc3RvcmFnZSBpZiBtaXNzaW5nLCByZXR1cm4gaXRcclxuICAgICAgICB2YXIgc3RvcmFnZSA9IGZ1bmN0aW9uKG9iail7XHJcbiAgICAgICAgICBpZiAoIWhhc093bihvYmosIGdsb2JhbElEKSlcclxuICAgICAgICAgICAgZGVmaW5lUHJvcGVydHkob2JqLCBnbG9iYWxJRCwgeyB2YWx1ZTogbmV3IEhhc2ggfSk7XHJcbiAgICAgICAgICByZXR1cm4gb2JqW2dsb2JhbElEXTtcclxuICAgICAgICB9O1xyXG5cclxuICAgICAgICBkZWZpbmUoT2JqZWN0LCBnZXRPd25Qcm9wZXJ0eU5hbWVzKTtcclxuICAgICAgfSBlbHNlIHtcclxuXHJcbiAgICAgICAgdmFyIHRvU3RyaW5nVG9TdHJpbmcgPSBmdW5jdGlvbihzKXtcclxuICAgICAgICAgIGZ1bmN0aW9uIHRvU3RyaW5nKCl7IHJldHVybiBzIH1cclxuICAgICAgICAgIHJldHVybiB0b1N0cmluZ1t0b1N0cmluZ19dID0gdG9TdHJpbmc7XHJcbiAgICAgICAgfShPYmplY3RbcHJvdG90eXBlX11bdG9TdHJpbmdfXSsnJyk7XHJcblxyXG4gICAgICAgIC8vIHN0b3JlIHRoZSB2YWx1ZXMgb24gYSBjdXN0b20gdmFsdWVPZiBpbiBvcmRlciB0byBoaWRlIHRoZW0gYnV0IHN0b3JlIHRoZW0gbG9jYWxseVxyXG4gICAgICAgIHZhciBzdG9yYWdlID0gZnVuY3Rpb24ob2JqKXtcclxuICAgICAgICAgIGlmIChoYXNPd24ob2JqLCB0b1N0cmluZ18pICYmIGdsb2JhbElEIGluIG9ialt0b1N0cmluZ19dKVxyXG4gICAgICAgICAgICByZXR1cm4gb2JqW3RvU3RyaW5nX11bZ2xvYmFsSURdO1xyXG5cclxuICAgICAgICAgIGlmICghKHRvU3RyaW5nXyBpbiBvYmopKVxyXG4gICAgICAgICAgICB0aHJvdyBuZXcgRXJyb3IoXCJDYW4ndCBzdG9yZSB2YWx1ZXMgZm9yIFwiK29iaik7XHJcblxyXG4gICAgICAgICAgdmFyIG9sZFRvU3RyaW5nID0gb2JqW3RvU3RyaW5nX107XHJcbiAgICAgICAgICBmdW5jdGlvbiB0b1N0cmluZygpeyByZXR1cm4gb2xkVG9TdHJpbmcuY2FsbCh0aGlzKSB9XHJcbiAgICAgICAgICBvYmpbdG9TdHJpbmdfXSA9IHRvU3RyaW5nO1xyXG4gICAgICAgICAgdG9TdHJpbmdbdG9TdHJpbmdfXSA9IHRvU3RyaW5nVG9TdHJpbmc7XHJcbiAgICAgICAgICByZXR1cm4gdG9TdHJpbmdbZ2xvYmFsSURdID0ge307XHJcbiAgICAgICAgfTtcclxuICAgICAgfVxyXG5cclxuXHJcblxyXG4gICAgICAvLyBzaGltIGZvciBbW01hcERhdGFdXSBmcm9tIGVzNiBzcGVjLCBhbmQgcHVsbHMgZG91YmxlIGR1dHkgYXMgV2Vha01hcCBzdG9yYWdlXHJcbiAgICAgIGZ1bmN0aW9uIE1hcERhdGEobmFtZSl7XHJcbiAgICAgICAgdmFyIHB1aWQgPSBjcmVhdGVVSUQoKSxcclxuICAgICAgICAgICAgaXVpZCA9IGNyZWF0ZVVJRCgpLFxyXG4gICAgICAgICAgICBzZWNyZXQgPSB7IHZhbHVlOiB1bmRlZmluZWQgfTtcclxuXHJcbiAgICAgICAgdmFyIGF0dGFjaCA9IGZ1bmN0aW9uKG9iail7XHJcbiAgICAgICAgICB2YXIgc3RvcmUgPSBzdG9yYWdlKG9iaik7XHJcbiAgICAgICAgICBpZiAoaGFzT3duKHN0b3JlLCBwdWlkKSlcclxuICAgICAgICAgICAgcmV0dXJuIHN0b3JlW3B1aWRdKHNlY3JldCk7XHJcblxyXG4gICAgICAgICAgdmFyIGxvY2tib3ggPSBuZXcgSGFzaDtcclxuICAgICAgICAgIGRlZmluZVByb3BlcnR5KGxvY2tib3gsIGl1aWQsIHNlY3JldCk7XHJcbiAgICAgICAgICBkZWZpbmVQcm9wZXJ0eShzdG9yZSwgcHVpZCwge1xyXG4gICAgICAgICAgICB2YWx1ZTogbmV3IEZ1bmN0aW9uKCdzJywgJ2wnLCBsb2NrZXIpKHNlY3JldCwgbG9ja2JveClcclxuICAgICAgICAgIH0pO1xyXG4gICAgICAgICAgcmV0dXJuIGxvY2tib3g7XHJcbiAgICAgICAgfTtcclxuXHJcbiAgICAgICAgdGhpcy5nZXQgPSBmdW5jdGlvbihvKXtcclxuICAgICAgICAgIHJldHVybiBhdHRhY2gobylbaXVpZF07XHJcbiAgICAgICAgfTtcclxuICAgICAgICB0aGlzLnNldCA9IGZ1bmN0aW9uKG8sIHYpe1xyXG4gICAgICAgICAgYXR0YWNoKG8pW2l1aWRdID0gdjtcclxuICAgICAgICB9O1xyXG5cclxuICAgICAgICBpZiAobmFtZSkge1xyXG4gICAgICAgICAgdGhpcy53cmFwID0gZnVuY3Rpb24obywgdil7XHJcbiAgICAgICAgICAgIHZhciBsb2NrYm94ID0gYXR0YWNoKG8pO1xyXG4gICAgICAgICAgICBpZiAobG9ja2JveFtpdWlkXSlcclxuICAgICAgICAgICAgICB0aHJvdyBuZXcgVHlwZUVycm9yKFwiT2JqZWN0IGlzIGFscmVhZHkgYSBcIiArIG5hbWUpO1xyXG4gICAgICAgICAgICBsb2NrYm94W2l1aWRdID0gdjtcclxuICAgICAgICAgIH07XHJcbiAgICAgICAgICB0aGlzLnVud3JhcCA9IGZ1bmN0aW9uKG8pe1xyXG4gICAgICAgICAgICB2YXIgc3RvcmFnZSA9IGF0dGFjaChvKVtpdWlkXTtcclxuICAgICAgICAgICAgaWYgKCFzdG9yYWdlKVxyXG4gICAgICAgICAgICAgIHRocm93IG5ldyBUeXBlRXJyb3IobmFtZSArIFwiIGlzIG5vdCBnZW5lcmljXCIpO1xyXG4gICAgICAgICAgICByZXR1cm4gc3RvcmFnZTtcclxuICAgICAgICAgIH07XHJcbiAgICAgICAgfVxyXG4gICAgICB9XHJcblxyXG4gICAgICByZXR1cm4gTWFwRGF0YTtcclxuICAgIH0oKSk7XHJcblxyXG4gIHZhciBleHBvcnRlciA9IChmdW5jdGlvbigpe1xyXG4gICAgLy8gW25hdGl2ZSBjb2RlXSBsb29rcyBzbGlnaHRseSBkaWZmZXJlbnQgaW4gZWFjaCBlbmdpbmVcclxuICAgIHZhciBzcmMgPSAoJycrT2JqZWN0KS5zcGxpdCgnT2JqZWN0Jyk7XHJcblxyXG4gICAgLy8gZmFrZSBbbmF0aXZlIGNvZGVdXHJcbiAgICBmdW5jdGlvbiB0b1N0cmluZygpe1xyXG4gICAgICByZXR1cm4gc3JjWzBdICsgbmFtZSh0aGlzKSArIHNyY1sxXTtcclxuICAgIH1cclxuXHJcbiAgICBkZWZpbmUodG9TdHJpbmcsIHRvU3RyaW5nKTtcclxuXHJcbiAgICAvLyBhdHRlbXB0IHRvIHVzZSBfX3Byb3RvX18gc28gdGhlIG1ldGhvZHMgZG9uJ3QgYWxsIGhhdmUgYW4gb3duIHRvU3RyaW5nXHJcbiAgICB2YXIgcHJlcEZ1bmN0aW9uID0geyBfX3Byb3RvX186IFtdIH0gaW5zdGFuY2VvZiBBcnJheVxyXG4gICAgICA/IGZ1bmN0aW9uKGZ1bmMpeyBmdW5jLl9fcHJvdG9fXyA9IHRvU3RyaW5nIH1cclxuICAgICAgOiBmdW5jdGlvbihmdW5jKXsgZGVmaW5lKGZ1bmMsIHRvU3RyaW5nKSB9O1xyXG5cclxuICAgIC8vIGFzc2VtYmxlIGFuIGFycmF5IG9mIGZ1bmN0aW9ucyBpbnRvIGEgZnVsbHkgZm9ybWVkIGNsYXNzXHJcbiAgICB2YXIgcHJlcGFyZSA9IGZ1bmN0aW9uKG1ldGhvZHMpe1xyXG4gICAgICB2YXIgQ3RvciA9IG1ldGhvZHMuc2hpZnQoKSxcclxuICAgICAgICAgIGJyYW5kID0gJ1tvYmplY3QgJyArIG5hbWUoQ3RvcikgKyAnXSc7XHJcblxyXG4gICAgICBmdW5jdGlvbiB0b1N0cmluZygpeyByZXR1cm4gYnJhbmQgfVxyXG4gICAgICBtZXRob2RzLnB1c2godG9TdHJpbmcpO1xyXG4gICAgICBwcmVwRnVuY3Rpb24oQ3Rvcik7XHJcblxyXG4gICAgICBmb3IgKHZhciBpPTA7IGkgPCBtZXRob2RzLmxlbmd0aDsgaSsrKSB7XHJcbiAgICAgICAgcHJlcEZ1bmN0aW9uKG1ldGhvZHNbaV0pO1xyXG4gICAgICAgIGRlZmluZShDdG9yW3Byb3RvdHlwZV9dLCBtZXRob2RzW2ldKTtcclxuICAgICAgfVxyXG5cclxuICAgICAgcmV0dXJuIEN0b3I7XHJcbiAgICB9O1xyXG5cclxuICAgIHJldHVybiBmdW5jdGlvbihuYW1lLCBpbml0KXtcclxuICAgICAgaWYgKG5hbWUgaW4gZXhwb3J0cylcclxuICAgICAgICByZXR1cm4gZXhwb3J0c1tuYW1lXTtcclxuXHJcbiAgICAgIHZhciBkYXRhID0gbmV3IE1hcERhdGEobmFtZSk7XHJcblxyXG4gICAgICByZXR1cm4gZXhwb3J0c1tuYW1lXSA9IHByZXBhcmUoaW5pdChcclxuICAgICAgICBmdW5jdGlvbihjb2xsZWN0aW9uLCB2YWx1ZSl7XHJcbiAgICAgICAgICBkYXRhLndyYXAoY29sbGVjdGlvbiwgdmFsdWUpO1xyXG4gICAgICAgIH0sXHJcbiAgICAgICAgZnVuY3Rpb24oY29sbGVjdGlvbil7XHJcbiAgICAgICAgICByZXR1cm4gZGF0YS51bndyYXAoY29sbGVjdGlvbik7XHJcbiAgICAgICAgfVxyXG4gICAgICApKTtcclxuICAgIH07XHJcbiAgfSgpKTtcclxuXHJcblxyXG4gIC8vIGluaXRpYWxpemUgY29sbGVjdGlvbiB3aXRoIGFuIGl0ZXJhYmxlLCBjdXJyZW50bHkgb25seSBzdXBwb3J0cyBmb3JFYWNoIGZ1bmN0aW9uXHJcbiAgdmFyIGluaXRpYWxpemUgPSBmdW5jdGlvbihpdGVyYWJsZSwgY2FsbGJhY2spe1xyXG4gICAgaWYgKGl0ZXJhYmxlICE9PSBudWxsICYmIHR5cGVvZiBpdGVyYWJsZSA9PT0gb2JqZWN0XyAmJiB0eXBlb2YgaXRlcmFibGUuZm9yRWFjaCA9PT0gZnVuY3Rpb25fKSB7XHJcbiAgICAgIGl0ZXJhYmxlLmZvckVhY2goZnVuY3Rpb24oaXRlbSwgaSl7XHJcbiAgICAgICAgaWYgKGlzQXJyYXkoaXRlbSkgJiYgaXRlbS5sZW5ndGggPT09IDIpXHJcbiAgICAgICAgICBjYWxsYmFjayhpdGVyYWJsZVtpXVswXSwgaXRlcmFibGVbaV1bMV0pO1xyXG4gICAgICAgIGVsc2VcclxuICAgICAgICAgIGNhbGxiYWNrKGl0ZXJhYmxlW2ldLCBpKTtcclxuICAgICAgfSk7XHJcbiAgICB9XHJcbiAgfVxyXG5cclxuICAvLyBhdHRlbXB0IHRvIGZpeCB0aGUgbmFtZSBvZiBcImRlbGV0ZV9cIiBtZXRob2RzLCBzaG91bGQgd29yayBpbiBWOCBhbmQgc3BpZGVybW9ua2V5XHJcbiAgdmFyIGZpeERlbGV0ZSA9IGZ1bmN0aW9uKGZ1bmMsIHNjb3BlTmFtZXMsIHNjb3BlVmFsdWVzKXtcclxuICAgIHRyeSB7XHJcbiAgICAgIHNjb3BlTmFtZXNbc2NvcGVOYW1lcy5sZW5ndGhdID0gKCdyZXR1cm4gJytmdW5jKS5yZXBsYWNlKCdlXycsICdcXFxcdTAwNjUnKTtcclxuICAgICAgcmV0dXJuIEZ1bmN0aW9uLmFwcGx5KDAsIHNjb3BlTmFtZXMpLmFwcGx5KDAsIHNjb3BlVmFsdWVzKTtcclxuICAgIH0gY2F0Y2ggKGUpIHtcclxuICAgICAgcmV0dXJuIGZ1bmM7XHJcbiAgICB9XHJcbiAgfVxyXG5cclxuICB2YXIgV00sIEhNLCBNO1xyXG5cclxuICAvLyAjIyMjIyMjIyMjIyMjIyNcclxuICAvLyAjIyMgV2Vha01hcCAjIyNcclxuICAvLyAjIyMjIyMjIyMjIyMjIyNcclxuXHJcbiAgV00gPSBidWlsdGluV2Vha01hcCA/IChleHBvcnRzLldlYWtNYXAgPSBnbG9iYWwuV2Vha01hcCkgOiBleHBvcnRlcignV2Vha01hcCcsIGZ1bmN0aW9uKHdyYXAsIHVud3JhcCl7XHJcbiAgICB2YXIgcHJvdG90eXBlID0gV2Vha01hcFtwcm90b3R5cGVfXTtcclxuICAgIHZhciB2YWxpZGF0ZSA9IGZ1bmN0aW9uKGtleSl7XHJcbiAgICAgIGlmIChrZXkgPT0gbnVsbCB8fCB0eXBlb2Yga2V5ICE9PSBvYmplY3RfICYmIHR5cGVvZiBrZXkgIT09IGZ1bmN0aW9uXylcclxuICAgICAgICB0aHJvdyBuZXcgVHlwZUVycm9yKFwiSW52YWxpZCBXZWFrTWFwIGtleVwiKTtcclxuICAgIH07XHJcblxyXG4gICAgLyoqXHJcbiAgICAgKiBAY2xhc3MgICAgICAgIFdlYWtNYXBcclxuICAgICAqIEBkZXNjcmlwdGlvbiAgQ29sbGVjdGlvbiB1c2luZyBvYmplY3RzIHdpdGggdW5pcXVlIGlkZW50aXRpZXMgYXMga2V5cyB0aGF0IGRpc2FsbG93cyBlbnVtZXJhdGlvblxyXG4gICAgICogICAgICAgICAgICAgICBhbmQgYWxsb3dzIGZvciBiZXR0ZXIgZ2FyYmFnZSBjb2xsZWN0aW9uLlxyXG4gICAgICogQHBhcmFtICAgICAgICB7SXRlcmFibGV9IFtpdGVyYWJsZV0gIEFuIGl0ZW0gdG8gcG9wdWxhdGUgdGhlIGNvbGxlY3Rpb24gd2l0aC5cclxuICAgICAqL1xyXG4gICAgZnVuY3Rpb24gV2Vha01hcChpdGVyYWJsZSl7XHJcbiAgICAgIGlmICh0aGlzID09PSBnbG9iYWwgfHwgdGhpcyA9PSBudWxsIHx8IHRoaXMgPT09IHByb3RvdHlwZSlcclxuICAgICAgICByZXR1cm4gbmV3IFdlYWtNYXAoaXRlcmFibGUpO1xyXG5cclxuICAgICAgd3JhcCh0aGlzLCBuZXcgTWFwRGF0YSk7XHJcblxyXG4gICAgICB2YXIgc2VsZiA9IHRoaXM7XHJcbiAgICAgIGl0ZXJhYmxlICYmIGluaXRpYWxpemUoaXRlcmFibGUsIGZ1bmN0aW9uKHZhbHVlLCBrZXkpe1xyXG4gICAgICAgIGNhbGwoc2V0LCBzZWxmLCB2YWx1ZSwga2V5KTtcclxuICAgICAgfSk7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPGdldD5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgUmV0cmlldmUgdGhlIHZhbHVlIGluIHRoZSBjb2xsZWN0aW9uIHRoYXQgbWF0Y2hlcyBrZXlcclxuICAgICAqIEBwYXJhbSAgICAgICAge0FueX0ga2V5XHJcbiAgICAgKiBAcmV0dXJuICAgICAgIHtBbnl9XHJcbiAgICAgKi9cclxuICAgIGZ1bmN0aW9uIGdldChrZXkpe1xyXG4gICAgICB2YWxpZGF0ZShrZXkpO1xyXG4gICAgICB2YXIgdmFsdWUgPSB1bndyYXAodGhpcykuZ2V0KGtleSk7XHJcbiAgICAgIHJldHVybiB2YWx1ZSA9PT0gdW5kZWZpbmVkXyA/IHVuZGVmaW5lZCA6IHZhbHVlO1xyXG4gICAgfVxyXG4gICAgLyoqXHJcbiAgICAgKiBAbWV0aG9kICAgICAgIDxzZXQ+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIEFkZCBvciB1cGRhdGUgYSBwYWlyIGluIHRoZSBjb2xsZWN0aW9uLiBFbmZvcmNlcyB1bmlxdWVuZXNzIGJ5IG92ZXJ3cml0aW5nLlxyXG4gICAgICogQHBhcmFtICAgICAgICB7QW55fSBrZXlcclxuICAgICAqIEBwYXJhbSAgICAgICAge0FueX0gdmFsXHJcbiAgICAgKiovXHJcbiAgICBmdW5jdGlvbiBzZXQoa2V5LCB2YWx1ZSl7XHJcbiAgICAgIHZhbGlkYXRlKGtleSk7XHJcbiAgICAgIC8vIHN0b3JlIGEgdG9rZW4gZm9yIGV4cGxpY2l0IHVuZGVmaW5lZCBzbyB0aGF0IFwiaGFzXCIgd29ya3MgY29ycmVjdGx5XHJcbiAgICAgIHVud3JhcCh0aGlzKS5zZXQoa2V5LCB2YWx1ZSA9PT0gdW5kZWZpbmVkID8gdW5kZWZpbmVkXyA6IHZhbHVlKTtcclxuICAgIH1cclxuICAgIC8qXHJcbiAgICAgKiBAbWV0aG9kICAgICAgIDxoYXM+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIENoZWNrIGlmIGtleSBpcyBpbiB0aGUgY29sbGVjdGlvblxyXG4gICAgICogQHBhcmFtICAgICAgICB7QW55fSBrZXlcclxuICAgICAqIEByZXR1cm4gICAgICAge0Jvb2xlYW59XHJcbiAgICAgKiovXHJcbiAgICBmdW5jdGlvbiBoYXMoa2V5KXtcclxuICAgICAgdmFsaWRhdGUoa2V5KTtcclxuICAgICAgcmV0dXJuIHVud3JhcCh0aGlzKS5nZXQoa2V5KSAhPT0gdW5kZWZpbmVkO1xyXG4gICAgfVxyXG4gICAgLyoqXHJcbiAgICAgKiBAbWV0aG9kICAgICAgIDxkZWxldGU+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIFJlbW92ZSBrZXkgYW5kIG1hdGNoaW5nIHZhbHVlIGlmIGZvdW5kXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICogQHJldHVybiAgICAgICB7Qm9vbGVhbn0gdHJ1ZSBpZiBpdGVtIHdhcyBpbiBjb2xsZWN0aW9uXHJcbiAgICAgKi9cclxuICAgIGZ1bmN0aW9uIGRlbGV0ZV8oa2V5KXtcclxuICAgICAgdmFsaWRhdGUoa2V5KTtcclxuICAgICAgdmFyIGRhdGEgPSB1bndyYXAodGhpcyk7XHJcblxyXG4gICAgICBpZiAoZGF0YS5nZXQoa2V5KSA9PT0gdW5kZWZpbmVkKVxyXG4gICAgICAgIHJldHVybiBmYWxzZTtcclxuXHJcbiAgICAgIGRhdGEuc2V0KGtleSwgdW5kZWZpbmVkKTtcclxuICAgICAgcmV0dXJuIHRydWU7XHJcbiAgICB9XHJcblxyXG4gICAgZGVsZXRlXyA9IGZpeERlbGV0ZShkZWxldGVfLCBbJ3ZhbGlkYXRlJywgJ3Vud3JhcCddLCBbdmFsaWRhdGUsIHVud3JhcF0pO1xyXG4gICAgcmV0dXJuIFtXZWFrTWFwLCBnZXQsIHNldCwgaGFzLCBkZWxldGVfXTtcclxuICB9KTtcclxuXHJcblxyXG4gIC8vICMjIyMjIyMjIyMjIyMjI1xyXG4gIC8vICMjIyBIYXNoTWFwICMjI1xyXG4gIC8vICMjIyMjIyMjIyMjIyMjI1xyXG5cclxuICBITSA9IGV4cG9ydGVyKCdIYXNoTWFwJywgZnVuY3Rpb24od3JhcCwgdW53cmFwKXtcclxuICAgIC8vIHNlcGFyYXRlIG51bWJlcnMsIHN0cmluZ3MsIGFuZCBhdG9tcyB0byBjb21wZW5zYXRlIGZvciBrZXkgY29lcmNpb24gdG8gc3RyaW5nXHJcblxyXG4gICAgdmFyIHByb3RvdHlwZSA9IEhhc2hNYXBbcHJvdG90eXBlX10sXHJcbiAgICAgICAgU1RSSU5HID0gMCwgTlVNQkVSID0gMSwgT1RIRVIgPSAyLFxyXG4gICAgICAgIG90aGVycyA9IHsgJ3RydWUnOiB0cnVlLCAnZmFsc2UnOiBmYWxzZSwgJ251bGwnOiBudWxsLCAwOiAtMCB9O1xyXG5cclxuICAgIHZhciBwcm90byA9IE1hdGgucmFuZG9tKCkudG9TdHJpbmcoMzYpLnNsaWNlKDIpO1xyXG5cclxuICAgIHZhciBjb2VyY2UgPSBmdW5jdGlvbihrZXkpe1xyXG4gICAgICByZXR1cm4ga2V5ID09PSAnX19wcm90b19fJyA/IHByb3RvIDoga2V5O1xyXG4gICAgfTtcclxuXHJcbiAgICB2YXIgdW5jb2VyY2UgPSBmdW5jdGlvbih0eXBlLCBrZXkpe1xyXG4gICAgICBzd2l0Y2ggKHR5cGUpIHtcclxuICAgICAgICBjYXNlIFNUUklORzogcmV0dXJuIGtleSA9PT0gcHJvdG8gPyAnX19wcm90b19fJyA6IGtleTtcclxuICAgICAgICBjYXNlIE5VTUJFUjogcmV0dXJuICtrZXk7XHJcbiAgICAgICAgY2FzZSBPVEhFUjogcmV0dXJuIG90aGVyc1trZXldO1xyXG4gICAgICB9XHJcbiAgICB9XHJcblxyXG5cclxuICAgIHZhciB2YWxpZGF0ZSA9IGZ1bmN0aW9uKGtleSl7XHJcbiAgICAgIGlmIChrZXkgPT0gbnVsbCkgcmV0dXJuIE9USEVSO1xyXG4gICAgICBzd2l0Y2ggKHR5cGVvZiBrZXkpIHtcclxuICAgICAgICBjYXNlICdib29sZWFuJzogcmV0dXJuIE9USEVSO1xyXG4gICAgICAgIGNhc2Ugc3RyaW5nXzogcmV0dXJuIFNUUklORztcclxuICAgICAgICAvLyBuZWdhdGl2ZSB6ZXJvIGhhcyB0byBiZSBleHBsaWNpdGx5IGFjY291bnRlZCBmb3JcclxuICAgICAgICBjYXNlICdudW1iZXInOiByZXR1cm4ga2V5ID09PSAwICYmIEluZmluaXR5IC8ga2V5ID09PSAtSW5maW5pdHkgPyBPVEhFUiA6IE5VTUJFUjtcclxuICAgICAgICBkZWZhdWx0OiB0aHJvdyBuZXcgVHlwZUVycm9yKFwiSW52YWxpZCBIYXNoTWFwIGtleVwiKTtcclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuICAgIC8qKlxyXG4gICAgICogQGNsYXNzICAgICAgICAgIEhhc2hNYXBcclxuICAgICAqIEBkZXNjcmlwdGlvbiAgICBDb2xsZWN0aW9uIHRoYXQgb25seSBhbGxvd3MgcHJpbWl0aXZlcyB0byBiZSBrZXlzLlxyXG4gICAgICogQHBhcmFtICAgICAgICAgIHtJdGVyYWJsZX0gW2l0ZXJhYmxlXSAgQW4gaXRlbSB0byBwb3B1bGF0ZSB0aGUgY29sbGVjdGlvbiB3aXRoLlxyXG4gICAgICovXHJcbiAgICBmdW5jdGlvbiBIYXNoTWFwKGl0ZXJhYmxlKXtcclxuICAgICAgaWYgKHRoaXMgPT09IGdsb2JhbCB8fCB0aGlzID09IG51bGwgfHwgdGhpcyA9PT0gcHJvdG90eXBlKVxyXG4gICAgICAgIHJldHVybiBuZXcgSGFzaE1hcChpdGVyYWJsZSk7XHJcblxyXG4gICAgICB3cmFwKHRoaXMsIHtcclxuICAgICAgICBzaXplOiAwLFxyXG4gICAgICAgIDA6IG5ldyBIYXNoLFxyXG4gICAgICAgIDE6IG5ldyBIYXNoLFxyXG4gICAgICAgIDI6IG5ldyBIYXNoXHJcbiAgICAgIH0pO1xyXG5cclxuICAgICAgdmFyIHNlbGYgPSB0aGlzO1xyXG4gICAgICBpdGVyYWJsZSAmJiBpbml0aWFsaXplKGl0ZXJhYmxlLCBmdW5jdGlvbih2YWx1ZSwga2V5KXtcclxuICAgICAgICBjYWxsKHNldCwgc2VsZiwgdmFsdWUsIGtleSk7XHJcbiAgICAgIH0pO1xyXG4gICAgfVxyXG4gICAgLyoqXHJcbiAgICAgKiBAbWV0aG9kICAgICAgIDxnZXQ+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIFJldHJpZXZlIHRoZSB2YWx1ZSBpbiB0aGUgY29sbGVjdGlvbiB0aGF0IG1hdGNoZXMga2V5XHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICogQHJldHVybiAgICAgICB7QW55fVxyXG4gICAgICovXHJcbiAgICBmdW5jdGlvbiBnZXQoa2V5KXtcclxuICAgICAgcmV0dXJuIHVud3JhcCh0aGlzKVt2YWxpZGF0ZShrZXkpXVtjb2VyY2Uoa2V5KV07XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPHNldD5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgQWRkIG9yIHVwZGF0ZSBhIHBhaXIgaW4gdGhlIGNvbGxlY3Rpb24uIEVuZm9yY2VzIHVuaXF1ZW5lc3MgYnkgb3ZlcndyaXRpbmcuXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICogQHBhcmFtICAgICAgICB7QW55fSB2YWxcclxuICAgICAqKi9cclxuICAgIGZ1bmN0aW9uIHNldChrZXksIHZhbHVlKXtcclxuICAgICAgdmFyIGl0ZW1zID0gdW53cmFwKHRoaXMpLFxyXG4gICAgICAgICAgZGF0YSA9IGl0ZW1zW3ZhbGlkYXRlKGtleSldO1xyXG5cclxuICAgICAga2V5ID0gY29lcmNlKGtleSk7XHJcbiAgICAgIGtleSBpbiBkYXRhIHx8IGl0ZW1zLnNpemUrKztcclxuICAgICAgZGF0YVtrZXldID0gdmFsdWU7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPGhhcz5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgQ2hlY2sgaWYga2V5IGV4aXN0cyBpbiB0aGUgY29sbGVjdGlvbi5cclxuICAgICAqIEBwYXJhbSAgICAgICAge0FueX0ga2V5XHJcbiAgICAgKiBAcmV0dXJuICAgICAgIHtCb29sZWFufSBpcyBpbiBjb2xsZWN0aW9uXHJcbiAgICAgKiovXHJcbiAgICBmdW5jdGlvbiBoYXMoa2V5KXtcclxuICAgICAgcmV0dXJuIGNvZXJjZShrZXkpIGluIHVud3JhcCh0aGlzKVt2YWxpZGF0ZShrZXkpXTtcclxuICAgIH1cclxuICAgIC8qKlxyXG4gICAgICogQG1ldGhvZCAgICAgICA8ZGVsZXRlPlxyXG4gICAgICogQGRlc2NyaXB0aW9uICBSZW1vdmUga2V5IGFuZCBtYXRjaGluZyB2YWx1ZSBpZiBmb3VuZFxyXG4gICAgICogQHBhcmFtICAgICAgICB7QW55fSBrZXlcclxuICAgICAqIEByZXR1cm4gICAgICAge0Jvb2xlYW59IHRydWUgaWYgaXRlbSB3YXMgaW4gY29sbGVjdGlvblxyXG4gICAgICovXHJcbiAgICBmdW5jdGlvbiBkZWxldGVfKGtleSl7XHJcbiAgICAgIHZhciBpdGVtcyA9IHVud3JhcCh0aGlzKSxcclxuICAgICAgICAgIGRhdGEgPSBpdGVtc1t2YWxpZGF0ZShrZXkpXTtcclxuXHJcbiAgICAgIGtleSA9IGNvZXJjZShrZXkpO1xyXG4gICAgICBpZiAoa2V5IGluIGRhdGEpIHtcclxuICAgICAgICBkZWxldGUgZGF0YVtrZXldO1xyXG4gICAgICAgIGl0ZW1zLnNpemUtLTtcclxuICAgICAgICByZXR1cm4gdHJ1ZTtcclxuICAgICAgfVxyXG5cclxuICAgICAgcmV0dXJuIGZhbHNlO1xyXG4gICAgfVxyXG4gICAgLyoqXHJcbiAgICAgKiBAbWV0aG9kICAgICAgIDxzaXplPlxyXG4gICAgICogQGRlc2NyaXB0aW9uICBSZXRyaWV2ZSB0aGUgYW1vdW50IG9mIGl0ZW1zIGluIHRoZSBjb2xsZWN0aW9uXHJcbiAgICAgKiBAcmV0dXJuICAgICAgIHtOdW1iZXJ9XHJcbiAgICAgKi9cclxuICAgIGZ1bmN0aW9uIHNpemUoKXtcclxuICAgICAgcmV0dXJuIHVud3JhcCh0aGlzKS5zaXplO1xyXG4gICAgfVxyXG4gICAgLyoqXHJcbiAgICAgKiBAbWV0aG9kICAgICAgIDxmb3JFYWNoPlxyXG4gICAgICogQGRlc2NyaXB0aW9uICBMb29wIHRocm91Z2ggdGhlIGNvbGxlY3Rpb24gcmFpc2luZyBjYWxsYmFjayBmb3IgZWFjaFxyXG4gICAgICogQHBhcmFtICAgICAgICB7RnVuY3Rpb259IGNhbGxiYWNrICBgY2FsbGJhY2sodmFsdWUsIGtleSlgXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtPYmplY3R9ICAgY29udGV4dCAgICBUaGUgYHRoaXNgIGJpbmRpbmcgZm9yIGNhbGxiYWNrcywgZGVmYXVsdCBudWxsXHJcbiAgICAgKi9cclxuICAgIGZ1bmN0aW9uIGZvckVhY2goY2FsbGJhY2ssIGNvbnRleHQpe1xyXG4gICAgICB2YXIgZGF0YSA9IHVud3JhcCh0aGlzKTtcclxuICAgICAgY29udGV4dCA9IGNvbnRleHQgPT0gbnVsbCA/IGdsb2JhbCA6IGNvbnRleHQ7XHJcbiAgICAgIGZvciAodmFyIGk9MDsgaSA8IDM7IGkrKylcclxuICAgICAgICBmb3IgKHZhciBrZXkgaW4gZGF0YVtpXSlcclxuICAgICAgICAgIGNhbGwoY2FsbGJhY2ssIGNvbnRleHQsIGRhdGFbaV1ba2V5XSwgdW5jb2VyY2UoaSwga2V5KSwgdGhpcyk7XHJcbiAgICB9XHJcblxyXG4gICAgZGVsZXRlXyA9IGZpeERlbGV0ZShkZWxldGVfLCBbJ3ZhbGlkYXRlJywgJ3Vud3JhcCcsICdjb2VyY2UnXSwgW3ZhbGlkYXRlLCB1bndyYXAsIGNvZXJjZV0pO1xyXG4gICAgcmV0dXJuIFtIYXNoTWFwLCBnZXQsIHNldCwgaGFzLCBkZWxldGVfLCBzaXplLCBmb3JFYWNoXTtcclxuICB9KTtcclxuXHJcblxyXG4gIC8vICMjIyMjIyMjIyMjXHJcbiAgLy8gIyMjIE1hcCAjIyNcclxuICAvLyAjIyMjIyMjIyMjI1xyXG5cclxuICAvLyBpZiBhIGZ1bGx5IGltcGxlbWVudGVkIE1hcCBleGlzdHMgdGhlbiB1c2UgaXRcclxuICBpZiAoJ01hcCcgaW4gZ2xvYmFsICYmICdmb3JFYWNoJyBpbiBnbG9iYWwuTWFwLnByb3RvdHlwZSkge1xyXG4gICAgTSA9IGV4cG9ydHMuTWFwID0gZ2xvYmFsLk1hcDtcclxuICB9IGVsc2Uge1xyXG4gICAgTSA9IGV4cG9ydGVyKCdNYXAnLCBmdW5jdGlvbih3cmFwLCB1bndyYXApe1xyXG4gICAgICAvLyBhdHRlbXB0IHRvIHVzZSBhbiBleGlzdGluZyBwYXJ0aWFsbHkgaW1wbGVtZW50ZWQgTWFwXHJcbiAgICAgIHZhciBCdWlsdGluTWFwID0gZ2xvYmFsLk1hcCxcclxuICAgICAgICAgIHByb3RvdHlwZSA9IE1hcFtwcm90b3R5cGVfXSxcclxuICAgICAgICAgIHdtID0gV01bcHJvdG90eXBlX10sXHJcbiAgICAgICAgICBobSA9IChCdWlsdGluTWFwIHx8IEhNKVtwcm90b3R5cGVfXSxcclxuICAgICAgICAgIG1nZXQgICAgPSBbY2FsbGJpbmQoaG0uZ2V0KSwgY2FsbGJpbmQod20uZ2V0KV0sXHJcbiAgICAgICAgICBtc2V0ICAgID0gW2NhbGxiaW5kKGhtLnNldCksIGNhbGxiaW5kKHdtLnNldCldLFxyXG4gICAgICAgICAgbWhhcyAgICA9IFtjYWxsYmluZChobS5oYXMpLCBjYWxsYmluZCh3bS5oYXMpXSxcclxuICAgICAgICAgIG1kZWxldGUgPSBbY2FsbGJpbmQoaG1bJ2RlbGV0ZSddKSwgY2FsbGJpbmQod21bJ2RlbGV0ZSddKV07XHJcblxyXG4gICAgICB2YXIgdHlwZSA9IEJ1aWx0aW5NYXBcclxuICAgICAgICA/IGZ1bmN0aW9uKCl7IHJldHVybiAwIH1cclxuICAgICAgICA6IGZ1bmN0aW9uKG8peyByZXR1cm4gKyh0eXBlb2YgbyA9PT0gb2JqZWN0XyA/IG8gIT09IG51bGwgOiB0eXBlb2YgbyA9PT0gZnVuY3Rpb25fKSB9XHJcblxyXG4gICAgICAvLyBpZiB3ZSBoYXZlIGEgYnVpbHRpbiBNYXAgd2UgY2FuIGxldCBpdCBkbyBtb3N0IG9mIHRoZSBoZWF2eSBsaWZ0aW5nXHJcbiAgICAgIHZhciBpbml0ID0gQnVpbHRpbk1hcFxyXG4gICAgICAgID8gZnVuY3Rpb24oKXsgcmV0dXJuIHsgMDogbmV3IEJ1aWx0aW5NYXAgfSB9XHJcbiAgICAgICAgOiBmdW5jdGlvbigpeyByZXR1cm4geyAwOiBuZXcgSE0sIDE6IG5ldyBXTSB9IH07XHJcblxyXG4gICAgICAvKipcclxuICAgICAgICogQGNsYXNzICAgICAgICAgTWFwXHJcbiAgICAgICAqIEBkZXNjcmlwdGlvbiAgIENvbGxlY3Rpb24gdGhhdCBhbGxvd3MgYW55IGtpbmQgb2YgdmFsdWUgdG8gYmUgYSBrZXkuXHJcbiAgICAgICAqIEBwYXJhbSAgICAgICAgIHtJdGVyYWJsZX0gW2l0ZXJhYmxlXSAgQW4gaXRlbSB0byBwb3B1bGF0ZSB0aGUgY29sbGVjdGlvbiB3aXRoLlxyXG4gICAgICAgKi9cclxuICAgICAgZnVuY3Rpb24gTWFwKGl0ZXJhYmxlKXtcclxuICAgICAgICBpZiAodGhpcyA9PT0gZ2xvYmFsIHx8IHRoaXMgPT0gbnVsbCB8fCB0aGlzID09PSBwcm90b3R5cGUpXHJcbiAgICAgICAgICByZXR1cm4gbmV3IE1hcChpdGVyYWJsZSk7XHJcblxyXG4gICAgICAgIHZhciBkYXRhID0gaW5pdCgpO1xyXG4gICAgICAgIGRhdGEua2V5cyA9IFtdO1xyXG4gICAgICAgIGRhdGEudmFsdWVzID0gW107XHJcbiAgICAgICAgd3JhcCh0aGlzLCBkYXRhKTtcclxuXHJcbiAgICAgICAgdmFyIHNlbGYgPSB0aGlzO1xyXG4gICAgICAgIGl0ZXJhYmxlICYmIGluaXRpYWxpemUoaXRlcmFibGUsIGZ1bmN0aW9uKHZhbHVlLCBrZXkpe1xyXG4gICAgICAgICAgY2FsbChzZXQsIHNlbGYsIHZhbHVlLCBrZXkpO1xyXG4gICAgICAgIH0pO1xyXG4gICAgICB9XHJcbiAgICAgIC8qKlxyXG4gICAgICAgKiBAbWV0aG9kICAgICAgIDxnZXQ+XHJcbiAgICAgICAqIEBkZXNjcmlwdGlvbiAgUmV0cmlldmUgdGhlIHZhbHVlIGluIHRoZSBjb2xsZWN0aW9uIHRoYXQgbWF0Y2hlcyBrZXlcclxuICAgICAgICogQHBhcmFtICAgICAgICB7QW55fSBrZXlcclxuICAgICAgICogQHJldHVybiAgICAgICB7QW55fVxyXG4gICAgICAgKi9cclxuICAgICAgZnVuY3Rpb24gZ2V0KGtleSl7XHJcbiAgICAgICAgdmFyIGRhdGEgPSB1bndyYXAodGhpcyksXHJcbiAgICAgICAgICAgIHQgPSB0eXBlKGtleSk7XHJcbiAgICAgICAgcmV0dXJuIGRhdGEudmFsdWVzW21nZXRbdF0oZGF0YVt0XSwga2V5KV07XHJcbiAgICAgIH1cclxuICAgICAgLyoqXHJcbiAgICAgICAqIEBtZXRob2QgICAgICAgPHNldD5cclxuICAgICAgICogQGRlc2NyaXB0aW9uICBBZGQgb3IgdXBkYXRlIGEgcGFpciBpbiB0aGUgY29sbGVjdGlvbi4gRW5mb3JjZXMgdW5pcXVlbmVzcyBieSBvdmVyd3JpdGluZy5cclxuICAgICAgICogQHBhcmFtICAgICAgICB7QW55fSBrZXlcclxuICAgICAgICogQHBhcmFtICAgICAgICB7QW55fSB2YWxcclxuICAgICAgICoqL1xyXG4gICAgICBmdW5jdGlvbiBzZXQoa2V5LCB2YWx1ZSl7XHJcbiAgICAgICAgdmFyIGRhdGEgPSB1bndyYXAodGhpcyksXHJcbiAgICAgICAgICAgIHQgPSB0eXBlKGtleSksXHJcbiAgICAgICAgICAgIGluZGV4ID0gbWdldFt0XShkYXRhW3RdLCBrZXkpO1xyXG5cclxuICAgICAgICBpZiAoaW5kZXggPT09IHVuZGVmaW5lZCkge1xyXG4gICAgICAgICAgbXNldFt0XShkYXRhW3RdLCBrZXksIGRhdGEua2V5cy5sZW5ndGgpO1xyXG4gICAgICAgICAgcHVzaChkYXRhLmtleXMsIGtleSk7XHJcbiAgICAgICAgICBwdXNoKGRhdGEudmFsdWVzLCB2YWx1ZSk7XHJcbiAgICAgICAgfSBlbHNlIHtcclxuICAgICAgICAgIGRhdGEua2V5c1tpbmRleF0gPSBrZXk7XHJcbiAgICAgICAgICBkYXRhLnZhbHVlc1tpbmRleF0gPSB2YWx1ZTtcclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuICAgICAgLyoqXHJcbiAgICAgICAqIEBtZXRob2QgICAgICAgPGhhcz5cclxuICAgICAgICogQGRlc2NyaXB0aW9uICBDaGVjayBpZiBrZXkgZXhpc3RzIGluIHRoZSBjb2xsZWN0aW9uLlxyXG4gICAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICAgKiBAcmV0dXJuICAgICAgIHtCb29sZWFufSBpcyBpbiBjb2xsZWN0aW9uXHJcbiAgICAgICAqKi9cclxuICAgICAgZnVuY3Rpb24gaGFzKGtleSl7XHJcbiAgICAgICAgdmFyIHQgPSB0eXBlKGtleSk7XHJcbiAgICAgICAgcmV0dXJuIG1oYXNbdF0odW53cmFwKHRoaXMpW3RdLCBrZXkpO1xyXG4gICAgICB9XHJcbiAgICAgIC8qKlxyXG4gICAgICAgKiBAbWV0aG9kICAgICAgIDxkZWxldGU+XHJcbiAgICAgICAqIEBkZXNjcmlwdGlvbiAgUmVtb3ZlIGtleSBhbmQgbWF0Y2hpbmcgdmFsdWUgaWYgZm91bmRcclxuICAgICAgICogQHBhcmFtICAgICAgICB7QW55fSBrZXlcclxuICAgICAgICogQHJldHVybiAgICAgICB7Qm9vbGVhbn0gdHJ1ZSBpZiBpdGVtIHdhcyBpbiBjb2xsZWN0aW9uXHJcbiAgICAgICAqL1xyXG4gICAgICBmdW5jdGlvbiBkZWxldGVfKGtleSl7XHJcbiAgICAgICAgdmFyIGRhdGEgPSB1bndyYXAodGhpcyksXHJcbiAgICAgICAgICAgIHQgPSB0eXBlKGtleSksXHJcbiAgICAgICAgICAgIGluZGV4ID0gbWdldFt0XShkYXRhW3RdLCBrZXkpO1xyXG5cclxuICAgICAgICBpZiAoaW5kZXggPT09IHVuZGVmaW5lZClcclxuICAgICAgICAgIHJldHVybiBmYWxzZTtcclxuXHJcbiAgICAgICAgbWRlbGV0ZVt0XShkYXRhW3RdLCBrZXkpO1xyXG4gICAgICAgIHNwbGljZShkYXRhLmtleXMsIGluZGV4LCAxKTtcclxuICAgICAgICBzcGxpY2UoZGF0YS52YWx1ZXMsIGluZGV4LCAxKTtcclxuICAgICAgICByZXR1cm4gdHJ1ZTtcclxuICAgICAgfVxyXG4gICAgICAvKipcclxuICAgICAgICogQG1ldGhvZCAgICAgICA8c2l6ZT5cclxuICAgICAgICogQGRlc2NyaXB0aW9uICBSZXRyaWV2ZSB0aGUgYW1vdW50IG9mIGl0ZW1zIGluIHRoZSBjb2xsZWN0aW9uXHJcbiAgICAgICAqIEByZXR1cm4gICAgICAge051bWJlcn1cclxuICAgICAgICovXHJcbiAgICAgIGZ1bmN0aW9uIHNpemUoKXtcclxuICAgICAgICByZXR1cm4gdW53cmFwKHRoaXMpLmtleXMubGVuZ3RoO1xyXG4gICAgICB9XHJcbiAgICAgIC8qKlxyXG4gICAgICAgKiBAbWV0aG9kICAgICAgIDxmb3JFYWNoPlxyXG4gICAgICAgKiBAZGVzY3JpcHRpb24gIExvb3AgdGhyb3VnaCB0aGUgY29sbGVjdGlvbiByYWlzaW5nIGNhbGxiYWNrIGZvciBlYWNoXHJcbiAgICAgICAqIEBwYXJhbSAgICAgICAge0Z1bmN0aW9ufSBjYWxsYmFjayAgYGNhbGxiYWNrKHZhbHVlLCBrZXkpYFxyXG4gICAgICAgKiBAcGFyYW0gICAgICAgIHtPYmplY3R9ICAgY29udGV4dCAgICBUaGUgYHRoaXNgIGJpbmRpbmcgZm9yIGNhbGxiYWNrcywgZGVmYXVsdCBudWxsXHJcbiAgICAgICAqL1xyXG4gICAgICBmdW5jdGlvbiBmb3JFYWNoKGNhbGxiYWNrLCBjb250ZXh0KXtcclxuICAgICAgICB2YXIgZGF0YSA9IHVud3JhcCh0aGlzKSxcclxuICAgICAgICAgICAga2V5cyA9IGRhdGEua2V5cyxcclxuICAgICAgICAgICAgdmFsdWVzID0gZGF0YS52YWx1ZXM7XHJcblxyXG4gICAgICAgIGNvbnRleHQgPSBjb250ZXh0ID09IG51bGwgPyBnbG9iYWwgOiBjb250ZXh0O1xyXG5cclxuICAgICAgICBmb3IgKHZhciBpPTAsIGxlbj1rZXlzLmxlbmd0aDsgaSA8IGxlbjsgaSsrKVxyXG4gICAgICAgICAgY2FsbChjYWxsYmFjaywgY29udGV4dCwgdmFsdWVzW2ldLCBrZXlzW2ldLCB0aGlzKTtcclxuICAgICAgfVxyXG5cclxuICAgICAgZGVsZXRlXyA9IGZpeERlbGV0ZShkZWxldGVfLFxyXG4gICAgICAgIFsndHlwZScsICd1bndyYXAnLCAnY2FsbCcsICdzcGxpY2UnXSxcclxuICAgICAgICBbdHlwZSwgdW53cmFwLCBjYWxsLCBzcGxpY2VdXHJcbiAgICAgICk7XHJcbiAgICAgIHJldHVybiBbTWFwLCBnZXQsIHNldCwgaGFzLCBkZWxldGVfLCBzaXplLCBmb3JFYWNoXTtcclxuICAgIH0pO1xyXG4gIH1cclxuXHJcblxyXG4gIC8vICMjIyMjIyMjIyMjXHJcbiAgLy8gIyMjIFNldCAjIyNcclxuICAvLyAjIyMjIyMjIyMjI1xyXG5cclxuICBleHBvcnRlcignU2V0JywgZnVuY3Rpb24od3JhcCwgdW53cmFwKXtcclxuICAgIHZhciBwcm90b3R5cGUgPSBTZXRbcHJvdG90eXBlX10sXHJcbiAgICAgICAgbSA9IE1bcHJvdG90eXBlX10sXHJcbiAgICAgICAgbXNpemUgPSBjYWxsYmluZChtLnNpemUpLFxyXG4gICAgICAgIG1mb3JFYWNoID0gY2FsbGJpbmQobS5mb3JFYWNoKSxcclxuICAgICAgICBtZ2V0ID0gY2FsbGJpbmQobS5nZXQpLFxyXG4gICAgICAgIG1zZXQgPSBjYWxsYmluZChtLnNldCksXHJcbiAgICAgICAgbWhhcyA9IGNhbGxiaW5kKG0uaGFzKSxcclxuICAgICAgICBtZGVsZXRlID0gY2FsbGJpbmQobVsnZGVsZXRlJ10pO1xyXG5cclxuICAgIC8qKlxyXG4gICAgICogQGNsYXNzICAgICAgICBTZXRcclxuICAgICAqIEBkZXNjcmlwdGlvbiAgQ29sbGVjdGlvbiBvZiB2YWx1ZXMgdGhhdCBlbmZvcmNlcyB1bmlxdWVuZXNzLlxyXG4gICAgICogQHBhcmFtICAgICAgICB7SXRlcmFibGV9IFtpdGVyYWJsZV0gIEFuIGl0ZW0gdG8gcG9wdWxhdGUgdGhlIGNvbGxlY3Rpb24gd2l0aC5cclxuICAgICAqKi9cclxuICAgIGZ1bmN0aW9uIFNldChpdGVyYWJsZSl7XHJcbiAgICAgIGlmICh0aGlzID09PSBnbG9iYWwgfHwgdGhpcyA9PSBudWxsIHx8IHRoaXMgPT09IHByb3RvdHlwZSlcclxuICAgICAgICByZXR1cm4gbmV3IFNldChpdGVyYWJsZSk7XHJcblxyXG4gICAgICB3cmFwKHRoaXMsIG5ldyBNKTtcclxuXHJcbiAgICAgIHZhciBzZWxmID0gdGhpcztcclxuICAgICAgaXRlcmFibGUgJiYgaW5pdGlhbGl6ZShpdGVyYWJsZSwgZnVuY3Rpb24odmFsdWUsIGtleSl7XHJcbiAgICAgICAgY2FsbChhZGQsIHNlbGYsIGtleSk7XHJcbiAgICAgIH0pO1xyXG4gICAgfVxyXG4gICAgLyoqXHJcbiAgICAgKiBAbWV0aG9kICAgICAgIDxhZGQ+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIEluc2VydCB2YWx1ZSBpZiBub3QgZm91bmQsIGVuZm9yY2luZyB1bmlxdWVuZXNzLlxyXG4gICAgICogQHBhcmFtICAgICAgICB7QW55fSB2YWxcclxuICAgICAqL1xyXG4gICAgZnVuY3Rpb24gYWRkKGtleSl7XHJcbiAgICAgIG1zZXQodW53cmFwKHRoaXMpLCBrZXksIGtleSk7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPGhhcz5cclxuICAgICAqIEBkZXNjcmlwdGlvbiAgQ2hlY2sgaWYga2V5IGV4aXN0cyBpbiB0aGUgY29sbGVjdGlvbi5cclxuICAgICAqIEBwYXJhbSAgICAgICAge0FueX0ga2V5XHJcbiAgICAgKiBAcmV0dXJuICAgICAgIHtCb29sZWFufSBpcyBpbiBjb2xsZWN0aW9uXHJcbiAgICAgKiovXHJcbiAgICBmdW5jdGlvbiBoYXMoa2V5KXtcclxuICAgICAgcmV0dXJuIG1oYXModW53cmFwKHRoaXMpLCBrZXkpO1xyXG4gICAgfVxyXG4gICAgLyoqXHJcbiAgICAgKiBAbWV0aG9kICAgICAgIDxkZWxldGU+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIFJlbW92ZSBrZXkgYW5kIG1hdGNoaW5nIHZhbHVlIGlmIGZvdW5kXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtBbnl9IGtleVxyXG4gICAgICogQHJldHVybiAgICAgICB7Qm9vbGVhbn0gdHJ1ZSBpZiBpdGVtIHdhcyBpbiBjb2xsZWN0aW9uXHJcbiAgICAgKi9cclxuICAgIGZ1bmN0aW9uIGRlbGV0ZV8oa2V5KXtcclxuICAgICAgcmV0dXJuIG1kZWxldGUodW53cmFwKHRoaXMpLCBrZXkpO1xyXG4gICAgfVxyXG4gICAgLyoqXHJcbiAgICAgKiBAbWV0aG9kICAgICAgIDxzaXplPlxyXG4gICAgICogQGRlc2NyaXB0aW9uICBSZXRyaWV2ZSB0aGUgYW1vdW50IG9mIGl0ZW1zIGluIHRoZSBjb2xsZWN0aW9uXHJcbiAgICAgKiBAcmV0dXJuICAgICAgIHtOdW1iZXJ9XHJcbiAgICAgKi9cclxuICAgIGZ1bmN0aW9uIHNpemUoKXtcclxuICAgICAgcmV0dXJuIG1zaXplKHVud3JhcCh0aGlzKSk7XHJcbiAgICB9XHJcbiAgICAvKipcclxuICAgICAqIEBtZXRob2QgICAgICAgPGZvckVhY2g+XHJcbiAgICAgKiBAZGVzY3JpcHRpb24gIExvb3AgdGhyb3VnaCB0aGUgY29sbGVjdGlvbiByYWlzaW5nIGNhbGxiYWNrIGZvciBlYWNoLiBJbmRleCBpcyBzaW1wbHkgdGhlIGNvdW50ZXIgZm9yIHRoZSBjdXJyZW50IGl0ZXJhdGlvbi5cclxuICAgICAqIEBwYXJhbSAgICAgICAge0Z1bmN0aW9ufSBjYWxsYmFjayAgYGNhbGxiYWNrKHZhbHVlLCBpbmRleClgXHJcbiAgICAgKiBAcGFyYW0gICAgICAgIHtPYmplY3R9ICAgY29udGV4dCAgICBUaGUgYHRoaXNgIGJpbmRpbmcgZm9yIGNhbGxiYWNrcywgZGVmYXVsdCBudWxsXHJcbiAgICAgKi9cclxuICAgIGZ1bmN0aW9uIGZvckVhY2goY2FsbGJhY2ssIGNvbnRleHQpe1xyXG4gICAgICB2YXIgaW5kZXggPSAwLFxyXG4gICAgICAgICAgc2VsZiA9IHRoaXM7XHJcbiAgICAgIG1mb3JFYWNoKHVud3JhcCh0aGlzLCBmdW5jdGlvbihrZXkpe1xyXG4gICAgICAgIGNhbGwoY2FsbGJhY2ssIHRoaXMsIGtleSwgaW5kZXgrKywgc2VsZik7XHJcbiAgICAgIH0sIGNvbnRleHQpKTtcclxuICAgIH1cclxuXHJcbiAgICBkZWxldGVfID0gZml4RGVsZXRlKGRlbGV0ZV8sIFsnbWRlbGV0ZScsICd1bndyYXAnXSwgW21kZWxldGUsIHVud3JhcF0pO1xyXG4gICAgcmV0dXJuIFtTZXQsIGFkZCwgaGFzLCBkZWxldGVfLCBzaXplLCBmb3JFYWNoXTtcclxuICB9KTtcclxufSgnc3RyaW5nJywgJ29iamVjdCcsICdmdW5jdGlvbicsICdwcm90b3R5cGUnLCAndG9TdHJpbmcnLFxyXG4gIEFycmF5LCBPYmplY3QsIEZ1bmN0aW9uLCBGdW5jdGlvbi5wcm90b3R5cGUsICgwLCBldmFsKSgndGhpcycpLFxyXG4gIHR5cGVvZiBleHBvcnRzID09PSAndW5kZWZpbmVkJyA/IHRoaXMgOiBleHBvcnRzLCB7fSk7XHJcbiIsIlwidXNlIHN0cmljdFwiXG5cbnZhciBzaGVsbCA9IHJlcXVpcmUoXCJnYW1lLXNoZWxsXCIpKCk7XG52YXIgTWFwID0gcmVxdWlyZShcImhhcm1vbnktY29sbGVjdGlvbnNcIikuTWFwO1xudmFyIFNldCA9IHJlcXVpcmUoXCJoYXJtb255LWNvbGxlY3Rpb25zXCIpLlNldDtcbnZhciBkdCA9IHJlcXVpcmUoXCJkZWxhdW5heS10cmlhbmd1bGF0ZVwiKTtcblxudmFyIHNxcnQgPSBNYXRoLnNxcnQ7XG52YXIgcG93ID0gTWF0aC5wb3c7XG52YXIgZmxvb3IgPSBNYXRoLmZsb29yO1xudmFyIHJhbmRvbSA9IE1hdGgucmFuZG9tO1xuXG5cbnZhciBuYlBvaW50cyA9IDEwMDtcbnZhciBuYkFudHMgPSAyMDA7XG52YXIgbmJDaXRpZXMgPSAxMDtcblxuLy8gZ2VuZXJhdGUgcGxhY2VkIHBvaW50c1xudmFyIGNpdHlMaXN0ID0gW107XG53aGlsZShjaXR5TGlzdC5sZW5ndGggPCBuYkNpdGllcyl7XG4gICAgICB2YXIgcmFuZG9tbnVtYmVyPU1hdGguY2VpbChNYXRoLnJhbmRvbSgpKm5iUG9pbnRzKVxuICAgICAgdmFyIGZvdW5kPWZhbHNlO1xuICAgICAgZm9yKHZhciBpPTA7aTxjaXR5TGlzdC5sZW5ndGg7aSsrKXtcbiAgICAgICAgICAgIGlmKGNpdHlMaXN0W2ldPT1yYW5kb21udW1iZXIpe2ZvdW5kPXRydWU7YnJlYWt9XG4gICAgICB9XG4gICAgICBpZighZm91bmQpY2l0eUxpc3RbY2l0eUxpc3QubGVuZ3RoXT1yYW5kb21udW1iZXI7XG59XG52YXIgY2l0eVNldCA9IG5ldyBTZXQoY2l0eUxpc3QpO1xuXG5cbmZ1bmN0aW9uIHNpZ24oeCkgeyByZXR1cm4geCA/IHggPCAwID8gLTEgOiAxIDogMDsgfVxuXG4vL0luaXRpYWxpemUgdHJpYW5ndWxhdGlvblxudmFyIHBvaW50cyA9IG5ldyBBcnJheShuYlBvaW50cylcbiAgICBmb3IodmFyIGk9MDsgaTxuYlBvaW50czsgKytpKSB7XG4gICAgICAgIHBvaW50c1tpXSA9IHtpZCA6IGksIHg6cmFuZG9tKCksIHk6cmFuZG9tKCl9O1xuICAgIH1cbnZhciBjZWxscyA9IGR0KHBvaW50cy5tYXAoZnVuY3Rpb24ocCl7cmV0dXJuIFtwLngsIHAueV19KSlcblxuLy8gY3JlYXRlIGVkZ2VzXG52YXIgbmV4dEVkZ2VzID0gbmV3IE1hcCgpO1xudmFyIGVkZ2VzID0gW107XG52YXIgcGVybXV0YXRpb25zID0gW1swLDFdLCBbMSwwXSwgWzAsMl0sIFsyLDBdLCBbMSwyXSwgWzIsMV1dO1xudmFyIG5iRWRnZXMgPSAwO1xuY2VsbHMuZm9yRWFjaChmdW5jdGlvbihjZWxsKXtcbiAgZm9yKHZhciBpPTA7IGk8cGVybXV0YXRpb25zLmxlbmd0aDsgKytpKXtcbiAgICB2YXIgcyA9IHBlcm11dGF0aW9uc1tpXVswXTtcbiAgICB2YXIgZCA9IHBlcm11dGF0aW9uc1tpXVsxXTtcbiAgICB2YXIgcHMgPSBwb2ludHNbY2VsbFtzXV07XG4gICAgdmFyIHBkID0gcG9pbnRzW2NlbGxbZF1dO1xuICAgIHZhciBlZGdlID0ge2lkIDogbmJFZGdlcyxcbiAgICAgICAgICBzb3VyY2U6IGNlbGxbc10sIFxuICAgICAgICAgIGRlc3RpbmF0aW9uOiBjZWxsW2RdLCBcbiAgICAgICAgICBkaXN0YW5jZSA6IHNxcnQoIHBvdyhwcy54IC0gcGQueCwgMikgKyBwb3cocHMueSAtIHBkLnksIDIpICksXG4gICAgICAgICAgZGlyZWN0aW9uIDogTWF0aC5hdGFuKChwZC55LXBzLnkpLyhwZC54LXBzLngpKSxcbiAgICAgICAgICBvcmllbnRhdGlvbiA6IHNpZ24oKHBkLngtcHMueCkpLFxuICAgICAgICAgIHBoZXJvbW9uIDogMFxuICAgICAgICAgIH07XG4gICAgdmFyIG5leHRzO1xuICAgIGlmKG5leHRFZGdlcy5oYXMocHMuaWQpKXtcbiAgICAgICAgbmV4dHMgPSBuZXh0RWRnZXMuZ2V0KHBzLmlkKTtcbiAgICAgICAgbmV4dHMucHVzaChlZGdlKTtcbiAgICB9IGVsc2Uge1xuICAgICAgICBuZXh0cyA9IFtlZGdlXTtcbiAgICB9XG4gICAgbmV4dEVkZ2VzLnNldChwcy5pZCwgbmV4dHMpO1xuICAgIGVkZ2VzLnB1c2goZWRnZSk7XG4gICAgbmJFZGdlcysrO1xuICB9XG4gIFxufSlcblxuLy8gaW5pdGlhbGl6ZSBhbnRzXG52YXIgcG9wdWxhdGlvbiA9IG5ldyBBcnJheShuYkFudHMpO1xudmFyIGksajtcbmZvciAoaSA9IDA7IGkgPCBuYkFudHM7IGkrKykge1xuICAgIC8vIHRha2UgYSByYW5kb20gZWRnZVxuICAgIHZhciBlZGdlID0gZWRnZXNbTWF0aC5mbG9vcihlZGdlcy5sZW5ndGgqcmFuZG9tKCkpXTtcbiAgICB2YXIgeCA9IHBvaW50c1tlZGdlLnNvdXJjZV0ueCBcbiAgICB2YXIgeSA9IHBvaW50c1tlZGdlLnNvdXJjZV0ueVxuICAgIHBvcHVsYXRpb25baV0gPSBuZXcgQW50KHgsIHksIGVkZ2UpO1xufVxuXG52YXIgY2FudmFzLCBjb250ZXh0O1xuXG5zaGVsbC5vbihcImluaXRcIiwgZnVuY3Rpb24oKSB7XG4gICAgY2FudmFzID0gZG9jdW1lbnQuY3JlYXRlRWxlbWVudChcImNhbnZhc1wiKTtcbiAgICBjYW52YXMud2lkdGggPSBzaGVsbC53aWR0aDtcbiAgICBjYW52YXMuaGVpZ2h0ID0gc2hlbGwuaGVpZ2h0O1xuICAgIGNvbnRleHQgPSBjYW52YXMuZ2V0Q29udGV4dChcIjJkXCIpO1xuICAgIHNoZWxsLmVsZW1lbnQuYXBwZW5kQ2hpbGQoY2FudmFzKTtcbn0pXG5cbnNoZWxsLm9uKFwicmVzaXplXCIsIGZ1bmN0aW9uKHcsIGgpIHtcbiAgICBjYW52YXMud2lkdGggPSB3O1xuICAgIGNhbnZhcy5oZWlnaHQgPSBoO1xufSlcblxuc2hlbGwub24oXCJyZW5kZXJcIiwgZnVuY3Rpb24oKSB7XG4gICAgdmFyIHcgPSBjYW52YXMud2lkdGg7XG4gICAgdmFyIGggPSBjYW52YXMuaGVpZ2h0O1xuICAgIHZhciBtb3VzZSA9IFtzaGVsbC5tb3VzZVgvdywgc2hlbGwubW91c2VZL2hdO1xuICAgIGNvbnRleHQuc2V0VHJhbnNmb3JtKHcsIDAsIDAsIGgsIDAsIDApO1xuICAgIGNvbnRleHQuZmlsbFN0eWxlID0gXCIjZmZmXCI7XG4gICAgY29udGV4dC5maWxsUmVjdCgwLDAsdyxoKTtcblxuICAgIC8vIGVkZ2VzXG4gICAgY29udGV4dC5zdHJva2VTdHlsZSA9IFwiIzAwMFwiO1xuICAgIGZvcih2YXIgaT0wOyBpPGVkZ2VzLmxlbmd0aDsgKytpKSB7XG4gICAgICAgIHZhciBlZGdlID0gZWRnZXNbaV07XG4gICAgICAgIGlmIChlZGdlLnBoZXJvbW9uICE9IDApe1xuICAgICAgICAgICAgY29udGV4dC5saW5lV2lkdGggPSAwLjAwMSAqIGVkZ2UucGhlcm9tb247XG4gICAgICAgIH1lbHNlIHtcbiAgICAgICAgICAgIGNvbnRleHQubGluZVdpZHRoID0gMC4wMDAwMTtcbiAgICAgICAgfVxuICAgICAgICBjb250ZXh0LmJlZ2luUGF0aCgpO1xuICAgICAgICBjb250ZXh0Lm1vdmVUbyhwb2ludHNbZWRnZS5zb3VyY2VdLngsIHBvaW50c1tlZGdlLnNvdXJjZV0ueSk7XG4gICAgICAgIGNvbnRleHQubGluZVRvKHBvaW50c1tlZGdlLmRlc3RpbmF0aW9uXS54LCBwb2ludHNbZWRnZS5kZXN0aW5hdGlvbl0ueSk7XG4gICAgICAgIGNvbnRleHQuc3Ryb2tlKCk7XG4gICAgfVxuXG4gICAgLy8gdmVydGljZXNcbiAgICBmb3IodmFyIGk9MDsgaTxwb2ludHMubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgY29udGV4dC5iZWdpblBhdGgoKVxuICAgICAgICB2YXIgcG9pbnQgPSBwb2ludHNbaV07XG4gICAgICAgIGlmIChjaXR5U2V0Lmhhcyhwb2ludC5pZCkpIHtcbiAgICAgICAgICAgIGNvbnRleHQuZmlsbFN0eWxlID0gXCIjMDEwMURGXCI7XG4gICAgICAgICAgICBjb250ZXh0LmFyYyhwb2ludC54LCBwb2ludC55LCAwLjAwNiwgMCwgMipNYXRoLlBJKTtcbiAgICAgICAgfVxuICAgICAgICBlbHNlIHtcbiAgICAgICAgICAgIGNvbnRleHQuZmlsbFN0eWxlID0gXCIjMDAwXCI7XG4gICAgICAgICAgICBjb250ZXh0LmFyYyhwb2ludHNbaV0ueCwgcG9pbnRzW2ldLnksIDAuMDAzLCAwLCAyKk1hdGguUEkpO1xuICAgICAgICB9XG4gICAgICAgIGNvbnRleHQuY2xvc2VQYXRoKCk7XG4gICAgICAgIGNvbnRleHQuZmlsbCgpO1xuICAgIH1cblxuICAgIC8vIG1vdmUgYW50c1xuICAgIGZvciAoaSA9IDA7IGkgPCBuYkFudHM7IGkrKykge1xuICAgICAgICBwb3B1bGF0aW9uW2ldLnRyYW5zaXQoKTtcbiAgICB9XG5cbiAgICAvLyBwaGVyb21vbiBldmFwb3JhdGlvblxuICAgIGZvciAoaSA9IDA7IGkgPCBlZGdlcy5sZW5ndGg7IGkrKykge1xuICAgICAgICBpZihlZGdlc1tpXS5waGVyb21vbiA+IDApe1xuICAgICAgICAgICAgZWRnZXNbaV0ucGhlcm9tb24gLT0gMC4wMDE7XG4gICAgICAgIH1cbiAgICB9XG5cblxuICAgIGZvcih2YXIgaT0wOyBpPHBvcHVsYXRpb24ubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgY29udGV4dC5iZWdpblBhdGgoKVxuICAgICAgICB2YXIgeCA9IHBvcHVsYXRpb25baV0ucG9zWCAvLysgMC4wMSpyYW5kb20oKTtcbiAgICAgICAgdmFyIHkgPSBwb3B1bGF0aW9uW2ldLnBvc1kgLy8rIDAuMDEqcmFuZG9tKCk7XG4gICAgICAgIGlmIChwb3B1bGF0aW9uW2ldLnN0YXRlID09PSBcInBoZXJvbW9uaW5nXCIpe2NvbnRleHQuZmlsbFN0eWxlID0gXCIjRkYwMDAwXCJ9XG4gICAgICAgIGVsc2Uge2NvbnRleHQuZmlsbFN0eWxlID0gXCIjNjEwQjBCXCJ9XG4gICAgICAgIGNvbnRleHQuYXJjKHgsIHksIDAuMDAzLCAwLCAyKk1hdGguUEkpXG4gICAgICAgIGNvbnRleHQuY2xvc2VQYXRoKClcbiAgICAgICAgY29udGV4dC5maWxsKClcbiAgICB9XG4gIFxufSlcblxuZnVuY3Rpb24gQW50KHgsIHksIGVkZ2UpIHsgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuICAgIHRoaXMucG9zWCA9IHg7ICAgICAgICAgICAgICAgIFxuICAgIHRoaXMucG9zWSA9IHk7XG4gICAgdGhpcy5lZGdlID0gZWRnZTtcbiAgICB0aGlzLnN0ZXAgPSAwO1xuICAgIHRoaXMuc3RhdGUgPSBcImZvcmFnZVwiO1xuICAgIHRoaXMudHJhbnNpdCA9IHN0YXRlbWFjaGluZTsgXG4gICAgdGhpcy5tb3ZlID0gbW92ZTtcbiAgICB0aGlzLmVkZ2VzID0gW107XG4gICAgdGhpcy5sYXN0Q2l0eSA9IHVuZGVmaW5lZDtcbn1cbi8vIGZvcmFnZTogdGhlIGFudCB3YW5kZXJzIGFyb3VuZCB3aXRob3V0IGFueSBwaGVyb21vbiBkZXBvc2l0aW9uXG4vLyBvbmNlIGl0IGZpbmRzIGEgY2l0eSwgaXQgc3RhcnRzIHJlbWVtYmVyaW5nIHRoZSBub2RlcyBpdCBnb2VzIHRocm91Z2hcbi8vIHdoZW4gaXQgZmluZHMgYW5vdGhlciBjaXR5LCBpdCBjb21wdXRlcyB0aGUgcGF0aCBsZW5ndGggYW5kIGFkZHMgcGhlcm9tb25zIG9uZSBlYWNoIGVkZ2VzXG4vLyBwcm9wb3J0aW9ubmFseSB0byB0aGUgc2hvcnRlc3RuZXNzIG9mIHRoZSBwYXRoXG4vLyBpdCByZXNldHMgdGhlIGxpc3Qgb2Ygbm9kZXMgYW5kIGNvbnRpbnVlc1xuLy8gd2hpbGUgZm9yYWdpbmcgdGhlIGFudCBjaG9zZXMgdGhlIHBhdGggd2l0aCBhIHBoZXJvbW9uIHByZWZlcmVuY2VcblxuZnVuY3Rpb24gbW92ZSgpIHtcbiAgICB2YXIgZWRnZUNoYW5nZWQ7XG4gICAgdmFyIGNpdHlSZWFjaGVkID0gZmFsc2U7XG4gICAgLy8gb24gZWRnZVxuICAgIGlmICh0aGlzLnN0ZXAgPCB0aGlzLmVkZ2UuZGlzdGFuY2Upe1xuICAgICAgICB0aGlzLnBvc1ggKz0gMC4wMDUqTWF0aC5jb3ModGhpcy5lZGdlLmRpcmVjdGlvbikqdGhpcy5lZGdlLm9yaWVudGF0aW9uO1xuICAgICAgICB0aGlzLnBvc1kgKz0gMC4wMDUqTWF0aC5zaW4odGhpcy5lZGdlLmRpcmVjdGlvbikqdGhpcy5lZGdlLm9yaWVudGF0aW9uO1xuICAgICAgICB0aGlzLnN0ZXAgKz0gMC4wMDU7XG4gICAgICAgIGVkZ2VDaGFuZ2VkID0gZmFsc2U7XG4gICAgLy8gb24gdmVydGV4XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdGhpcy5zdGVwID0gMDtcbiAgICAgICAgdGhpcy5wb3NYID0gcG9pbnRzW3RoaXMuZWRnZS5kZXN0aW5hdGlvbl0ueDtcbiAgICAgICAgdGhpcy5wb3NZID0gcG9pbnRzW3RoaXMuZWRnZS5kZXN0aW5hdGlvbl0ueTtcbiAgICAgICAgdmFyIHBvc3NpYmxlRWRnZXMgPSBuZXh0RWRnZXMuZ2V0KHRoaXMuZWRnZS5kZXN0aW5hdGlvbik7XG4gICAgICAgIC8vIGZsaXAgYSBjb2luIGFuZCBlaXRoZXIgdGFrZSB0aGUgc21lbGxpZXN0IHBhdGggb2YgYSByYW5kb20gb25lXG4gICAgICAgIGlmIChyYW5kb20oKSA+IDAuNSl7XG4gICAgICAgICAgICB2YXIgc21lbGxzID0gcG9zc2libGVFZGdlcy5tYXAoZnVuY3Rpb24oZSl7cmV0dXJuIGUucGhlcm9tb259KTtcbiAgICAgICAgICAgIHZhciBpbmRleCA9IHNtZWxscy5pbmRleE9mKE1hdGgubWF4LmFwcGx5KE1hdGgsIHNtZWxscykpO1xuICAgICAgICAgICAgdGhpcy5lZGdlID0gcG9zc2libGVFZGdlc1tpbmRleF07XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB0aGlzLmVkZ2UgPSBwb3NzaWJsZUVkZ2VzW2Zsb29yKHJhbmRvbSgpKnBvc3NpYmxlRWRnZXMubGVuZ3RoKV07XG4gICAgICAgIH1cbiAgICAgICAgY2l0eVJlYWNoZWQgPSBjaXR5U2V0Lmhhcyh0aGlzLmVkZ2Uuc291cmNlKTtcbiAgICAgICAgZWRnZUNoYW5nZWQgPSB0cnVlO1xuICAgIH1cbiAgICByZXR1cm4ge2NpdHlSZWFjaGVkOiBjaXR5UmVhY2hlZCwgZWRnZUNoYW5nZWQ6IGVkZ2VDaGFuZ2VkfTtcbn1cblxuXG5mdW5jdGlvbiBzdGF0ZW1hY2hpbmUoKSB7XG4gICAgc3dpdGNoICh0aGlzLnN0YXRlKSB7XG4gICAgICAgIGNhc2UgXCJmb3JhZ2VcIjpcbiAgICAgICAgICAgIHZhciByZXMgPSB0aGlzLm1vdmUoKTtcbiAgICAgICAgICAgIGlmIChyZXMuY2l0eVJlYWNoZWQpIHtcbiAgICAgICAgICAgICAgICB0aGlzLnN0YXRlID0gXCJwaGVyb21vbmluZ1wiO1xuICAgICAgICAgICAgICAgIHRoaXMubGFzdENpdHkgPSB0aGlzLmVkZ2Uuc291cmNlO1xuICAgICAgICAgICAgfTtcbiAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICBjYXNlIFwicGhlcm9tb25pbmdcIjpcbiAgICAgICAgICAgIHZhciByZXMgPSB0aGlzLm1vdmUoKTtcbiAgICAgICAgICAgIGlmIChyZXMuZWRnZUNoYW5nZWQpIHtcbiAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzLnB1c2godGhpcy5lZGdlKTtcbiAgICAgICAgICAgICAgICAvLyBmb3VuZCBhIGNpdHlcbiAgICAgICAgICAgICAgICBpZiAocmVzLmNpdHlSZWFjaGVkICYmICh0aGlzLmVkZ2Uuc291cmNlICE9IHRoaXMubGFzdENpdHkpICl7XG4gICAgICAgICAgICAgICAgICAgIC8vIGNvbXB1dGUgdGhlIGxlbmd0aCBvZiB0aGUgcGF0aFxuICAgICAgICAgICAgICAgICAgICB2YXIgcGF0aExlbmd0aCA9IHRoaXMuZWRnZXMubWFwKGZ1bmN0aW9uKGUpe3JldHVybiBlLmRpc3RhbmNlfSkucmVkdWNlKGZ1bmN0aW9uKGEsYil7cmV0dXJuIGEgKyBifSk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBkZWx0YVBoZXJvbW9uZSA9IDEvcGF0aExlbmd0aDtcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5lZGdlcy5mb3JFYWNoKGZ1bmN0aW9uKGUpe2UucGhlcm9tb24gKz0gZGVsdGFQaGVyb21vbmV9KTtcbiAgICAgICAgICAgICAgICAgICAgLy8gY29uc29sZS5sb2coZGVsdGFQaGVyb21vbmUsIHRoaXMuZWRnZXMpO1xuICAgICAgICAgICAgICAgICAgICB0aGlzLmVkZ2VzID0gW3RoaXMuZWRnZV07XG4gICAgICAgICAgICAgICAgICAgIHRoaXMubGFzdENpdHkgPSB0aGlzLmVkZ2Uuc291cmNlO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgICBicmVhaztcblxuICAgIH1cbn0iLCIvLyBDb3B5cmlnaHQgSm95ZW50LCBJbmMuIGFuZCBvdGhlciBOb2RlIGNvbnRyaWJ1dG9ycy5cbi8vXG4vLyBQZXJtaXNzaW9uIGlzIGhlcmVieSBncmFudGVkLCBmcmVlIG9mIGNoYXJnZSwgdG8gYW55IHBlcnNvbiBvYnRhaW5pbmcgYVxuLy8gY29weSBvZiB0aGlzIHNvZnR3YXJlIGFuZCBhc3NvY2lhdGVkIGRvY3VtZW50YXRpb24gZmlsZXMgKHRoZVxuLy8gXCJTb2Z0d2FyZVwiKSwgdG8gZGVhbCBpbiB0aGUgU29mdHdhcmUgd2l0aG91dCByZXN0cmljdGlvbiwgaW5jbHVkaW5nXG4vLyB3aXRob3V0IGxpbWl0YXRpb24gdGhlIHJpZ2h0cyB0byB1c2UsIGNvcHksIG1vZGlmeSwgbWVyZ2UsIHB1Ymxpc2gsXG4vLyBkaXN0cmlidXRlLCBzdWJsaWNlbnNlLCBhbmQvb3Igc2VsbCBjb3BpZXMgb2YgdGhlIFNvZnR3YXJlLCBhbmQgdG8gcGVybWl0XG4vLyBwZXJzb25zIHRvIHdob20gdGhlIFNvZnR3YXJlIGlzIGZ1cm5pc2hlZCB0byBkbyBzbywgc3ViamVjdCB0byB0aGVcbi8vIGZvbGxvd2luZyBjb25kaXRpb25zOlxuLy9cbi8vIFRoZSBhYm92ZSBjb3B5cmlnaHQgbm90aWNlIGFuZCB0aGlzIHBlcm1pc3Npb24gbm90aWNlIHNoYWxsIGJlIGluY2x1ZGVkXG4vLyBpbiBhbGwgY29waWVzIG9yIHN1YnN0YW50aWFsIHBvcnRpb25zIG9mIHRoZSBTb2Z0d2FyZS5cbi8vXG4vLyBUSEUgU09GVFdBUkUgSVMgUFJPVklERUQgXCJBUyBJU1wiLCBXSVRIT1VUIFdBUlJBTlRZIE9GIEFOWSBLSU5ELCBFWFBSRVNTXG4vLyBPUiBJTVBMSUVELCBJTkNMVURJTkcgQlVUIE5PVCBMSU1JVEVEIFRPIFRIRSBXQVJSQU5USUVTIE9GXG4vLyBNRVJDSEFOVEFCSUxJVFksIEZJVE5FU1MgRk9SIEEgUEFSVElDVUxBUiBQVVJQT1NFIEFORCBOT05JTkZSSU5HRU1FTlQuIElOXG4vLyBOTyBFVkVOVCBTSEFMTCBUSEUgQVVUSE9SUyBPUiBDT1BZUklHSFQgSE9MREVSUyBCRSBMSUFCTEUgRk9SIEFOWSBDTEFJTSxcbi8vIERBTUFHRVMgT1IgT1RIRVIgTElBQklMSVRZLCBXSEVUSEVSIElOIEFOIEFDVElPTiBPRiBDT05UUkFDVCwgVE9SVCBPUlxuLy8gT1RIRVJXSVNFLCBBUklTSU5HIEZST00sIE9VVCBPRiBPUiBJTiBDT05ORUNUSU9OIFdJVEggVEhFIFNPRlRXQVJFIE9SIFRIRVxuLy8gVVNFIE9SIE9USEVSIERFQUxJTkdTIElOIFRIRSBTT0ZUV0FSRS5cblxuZnVuY3Rpb24gRXZlbnRFbWl0dGVyKCkge1xuICB0aGlzLl9ldmVudHMgPSB0aGlzLl9ldmVudHMgfHwge307XG4gIHRoaXMuX21heExpc3RlbmVycyA9IHRoaXMuX21heExpc3RlbmVycyB8fCB1bmRlZmluZWQ7XG59XG5tb2R1bGUuZXhwb3J0cyA9IEV2ZW50RW1pdHRlcjtcblxuLy8gQmFja3dhcmRzLWNvbXBhdCB3aXRoIG5vZGUgMC4xMC54XG5FdmVudEVtaXR0ZXIuRXZlbnRFbWl0dGVyID0gRXZlbnRFbWl0dGVyO1xuXG5FdmVudEVtaXR0ZXIucHJvdG90eXBlLl9ldmVudHMgPSB1bmRlZmluZWQ7XG5FdmVudEVtaXR0ZXIucHJvdG90eXBlLl9tYXhMaXN0ZW5lcnMgPSB1bmRlZmluZWQ7XG5cbi8vIEJ5IGRlZmF1bHQgRXZlbnRFbWl0dGVycyB3aWxsIHByaW50IGEgd2FybmluZyBpZiBtb3JlIHRoYW4gMTAgbGlzdGVuZXJzIGFyZVxuLy8gYWRkZWQgdG8gaXQuIFRoaXMgaXMgYSB1c2VmdWwgZGVmYXVsdCB3aGljaCBoZWxwcyBmaW5kaW5nIG1lbW9yeSBsZWFrcy5cbkV2ZW50RW1pdHRlci5kZWZhdWx0TWF4TGlzdGVuZXJzID0gMTA7XG5cbi8vIE9idmlvdXNseSBub3QgYWxsIEVtaXR0ZXJzIHNob3VsZCBiZSBsaW1pdGVkIHRvIDEwLiBUaGlzIGZ1bmN0aW9uIGFsbG93c1xuLy8gdGhhdCB0byBiZSBpbmNyZWFzZWQuIFNldCB0byB6ZXJvIGZvciB1bmxpbWl0ZWQuXG5FdmVudEVtaXR0ZXIucHJvdG90eXBlLnNldE1heExpc3RlbmVycyA9IGZ1bmN0aW9uKG4pIHtcbiAgaWYgKCFpc051bWJlcihuKSB8fCBuIDwgMCB8fCBpc05hTihuKSlcbiAgICB0aHJvdyBUeXBlRXJyb3IoJ24gbXVzdCBiZSBhIHBvc2l0aXZlIG51bWJlcicpO1xuICB0aGlzLl9tYXhMaXN0ZW5lcnMgPSBuO1xuICByZXR1cm4gdGhpcztcbn07XG5cbkV2ZW50RW1pdHRlci5wcm90b3R5cGUuZW1pdCA9IGZ1bmN0aW9uKHR5cGUpIHtcbiAgdmFyIGVyLCBoYW5kbGVyLCBsZW4sIGFyZ3MsIGksIGxpc3RlbmVycztcblxuICBpZiAoIXRoaXMuX2V2ZW50cylcbiAgICB0aGlzLl9ldmVudHMgPSB7fTtcblxuICAvLyBJZiB0aGVyZSBpcyBubyAnZXJyb3InIGV2ZW50IGxpc3RlbmVyIHRoZW4gdGhyb3cuXG4gIGlmICh0eXBlID09PSAnZXJyb3InKSB7XG4gICAgaWYgKCF0aGlzLl9ldmVudHMuZXJyb3IgfHxcbiAgICAgICAgKGlzT2JqZWN0KHRoaXMuX2V2ZW50cy5lcnJvcikgJiYgIXRoaXMuX2V2ZW50cy5lcnJvci5sZW5ndGgpKSB7XG4gICAgICBlciA9IGFyZ3VtZW50c1sxXTtcbiAgICAgIGlmIChlciBpbnN0YW5jZW9mIEVycm9yKSB7XG4gICAgICAgIHRocm93IGVyOyAvLyBVbmhhbmRsZWQgJ2Vycm9yJyBldmVudFxuICAgICAgfSBlbHNlIHtcbiAgICAgICAgdGhyb3cgVHlwZUVycm9yKCdVbmNhdWdodCwgdW5zcGVjaWZpZWQgXCJlcnJvclwiIGV2ZW50LicpO1xuICAgICAgfVxuICAgICAgcmV0dXJuIGZhbHNlO1xuICAgIH1cbiAgfVxuXG4gIGhhbmRsZXIgPSB0aGlzLl9ldmVudHNbdHlwZV07XG5cbiAgaWYgKGlzVW5kZWZpbmVkKGhhbmRsZXIpKVxuICAgIHJldHVybiBmYWxzZTtcblxuICBpZiAoaXNGdW5jdGlvbihoYW5kbGVyKSkge1xuICAgIHN3aXRjaCAoYXJndW1lbnRzLmxlbmd0aCkge1xuICAgICAgLy8gZmFzdCBjYXNlc1xuICAgICAgY2FzZSAxOlxuICAgICAgICBoYW5kbGVyLmNhbGwodGhpcyk7XG4gICAgICAgIGJyZWFrO1xuICAgICAgY2FzZSAyOlxuICAgICAgICBoYW5kbGVyLmNhbGwodGhpcywgYXJndW1lbnRzWzFdKTtcbiAgICAgICAgYnJlYWs7XG4gICAgICBjYXNlIDM6XG4gICAgICAgIGhhbmRsZXIuY2FsbCh0aGlzLCBhcmd1bWVudHNbMV0sIGFyZ3VtZW50c1syXSk7XG4gICAgICAgIGJyZWFrO1xuICAgICAgLy8gc2xvd2VyXG4gICAgICBkZWZhdWx0OlxuICAgICAgICBsZW4gPSBhcmd1bWVudHMubGVuZ3RoO1xuICAgICAgICBhcmdzID0gbmV3IEFycmF5KGxlbiAtIDEpO1xuICAgICAgICBmb3IgKGkgPSAxOyBpIDwgbGVuOyBpKyspXG4gICAgICAgICAgYXJnc1tpIC0gMV0gPSBhcmd1bWVudHNbaV07XG4gICAgICAgIGhhbmRsZXIuYXBwbHkodGhpcywgYXJncyk7XG4gICAgfVxuICB9IGVsc2UgaWYgKGlzT2JqZWN0KGhhbmRsZXIpKSB7XG4gICAgbGVuID0gYXJndW1lbnRzLmxlbmd0aDtcbiAgICBhcmdzID0gbmV3IEFycmF5KGxlbiAtIDEpO1xuICAgIGZvciAoaSA9IDE7IGkgPCBsZW47IGkrKylcbiAgICAgIGFyZ3NbaSAtIDFdID0gYXJndW1lbnRzW2ldO1xuXG4gICAgbGlzdGVuZXJzID0gaGFuZGxlci5zbGljZSgpO1xuICAgIGxlbiA9IGxpc3RlbmVycy5sZW5ndGg7XG4gICAgZm9yIChpID0gMDsgaSA8IGxlbjsgaSsrKVxuICAgICAgbGlzdGVuZXJzW2ldLmFwcGx5KHRoaXMsIGFyZ3MpO1xuICB9XG5cbiAgcmV0dXJuIHRydWU7XG59O1xuXG5FdmVudEVtaXR0ZXIucHJvdG90eXBlLmFkZExpc3RlbmVyID0gZnVuY3Rpb24odHlwZSwgbGlzdGVuZXIpIHtcbiAgdmFyIG07XG5cbiAgaWYgKCFpc0Z1bmN0aW9uKGxpc3RlbmVyKSlcbiAgICB0aHJvdyBUeXBlRXJyb3IoJ2xpc3RlbmVyIG11c3QgYmUgYSBmdW5jdGlvbicpO1xuXG4gIGlmICghdGhpcy5fZXZlbnRzKVxuICAgIHRoaXMuX2V2ZW50cyA9IHt9O1xuXG4gIC8vIFRvIGF2b2lkIHJlY3Vyc2lvbiBpbiB0aGUgY2FzZSB0aGF0IHR5cGUgPT09IFwibmV3TGlzdGVuZXJcIiEgQmVmb3JlXG4gIC8vIGFkZGluZyBpdCB0byB0aGUgbGlzdGVuZXJzLCBmaXJzdCBlbWl0IFwibmV3TGlzdGVuZXJcIi5cbiAgaWYgKHRoaXMuX2V2ZW50cy5uZXdMaXN0ZW5lcilcbiAgICB0aGlzLmVtaXQoJ25ld0xpc3RlbmVyJywgdHlwZSxcbiAgICAgICAgICAgICAgaXNGdW5jdGlvbihsaXN0ZW5lci5saXN0ZW5lcikgP1xuICAgICAgICAgICAgICBsaXN0ZW5lci5saXN0ZW5lciA6IGxpc3RlbmVyKTtcblxuICBpZiAoIXRoaXMuX2V2ZW50c1t0eXBlXSlcbiAgICAvLyBPcHRpbWl6ZSB0aGUgY2FzZSBvZiBvbmUgbGlzdGVuZXIuIERvbid0IG5lZWQgdGhlIGV4dHJhIGFycmF5IG9iamVjdC5cbiAgICB0aGlzLl9ldmVudHNbdHlwZV0gPSBsaXN0ZW5lcjtcbiAgZWxzZSBpZiAoaXNPYmplY3QodGhpcy5fZXZlbnRzW3R5cGVdKSlcbiAgICAvLyBJZiB3ZSd2ZSBhbHJlYWR5IGdvdCBhbiBhcnJheSwganVzdCBhcHBlbmQuXG4gICAgdGhpcy5fZXZlbnRzW3R5cGVdLnB1c2gobGlzdGVuZXIpO1xuICBlbHNlXG4gICAgLy8gQWRkaW5nIHRoZSBzZWNvbmQgZWxlbWVudCwgbmVlZCB0byBjaGFuZ2UgdG8gYXJyYXkuXG4gICAgdGhpcy5fZXZlbnRzW3R5cGVdID0gW3RoaXMuX2V2ZW50c1t0eXBlXSwgbGlzdGVuZXJdO1xuXG4gIC8vIENoZWNrIGZvciBsaXN0ZW5lciBsZWFrXG4gIGlmIChpc09iamVjdCh0aGlzLl9ldmVudHNbdHlwZV0pICYmICF0aGlzLl9ldmVudHNbdHlwZV0ud2FybmVkKSB7XG4gICAgdmFyIG07XG4gICAgaWYgKCFpc1VuZGVmaW5lZCh0aGlzLl9tYXhMaXN0ZW5lcnMpKSB7XG4gICAgICBtID0gdGhpcy5fbWF4TGlzdGVuZXJzO1xuICAgIH0gZWxzZSB7XG4gICAgICBtID0gRXZlbnRFbWl0dGVyLmRlZmF1bHRNYXhMaXN0ZW5lcnM7XG4gICAgfVxuXG4gICAgaWYgKG0gJiYgbSA+IDAgJiYgdGhpcy5fZXZlbnRzW3R5cGVdLmxlbmd0aCA+IG0pIHtcbiAgICAgIHRoaXMuX2V2ZW50c1t0eXBlXS53YXJuZWQgPSB0cnVlO1xuICAgICAgY29uc29sZS5lcnJvcignKG5vZGUpIHdhcm5pbmc6IHBvc3NpYmxlIEV2ZW50RW1pdHRlciBtZW1vcnkgJyArXG4gICAgICAgICAgICAgICAgICAgICdsZWFrIGRldGVjdGVkLiAlZCBsaXN0ZW5lcnMgYWRkZWQuICcgK1xuICAgICAgICAgICAgICAgICAgICAnVXNlIGVtaXR0ZXIuc2V0TWF4TGlzdGVuZXJzKCkgdG8gaW5jcmVhc2UgbGltaXQuJyxcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5fZXZlbnRzW3R5cGVdLmxlbmd0aCk7XG4gICAgICBpZiAodHlwZW9mIGNvbnNvbGUudHJhY2UgPT09ICdmdW5jdGlvbicpIHtcbiAgICAgICAgLy8gbm90IHN1cHBvcnRlZCBpbiBJRSAxMFxuICAgICAgICBjb25zb2xlLnRyYWNlKCk7XG4gICAgICB9XG4gICAgfVxuICB9XG5cbiAgcmV0dXJuIHRoaXM7XG59O1xuXG5FdmVudEVtaXR0ZXIucHJvdG90eXBlLm9uID0gRXZlbnRFbWl0dGVyLnByb3RvdHlwZS5hZGRMaXN0ZW5lcjtcblxuRXZlbnRFbWl0dGVyLnByb3RvdHlwZS5vbmNlID0gZnVuY3Rpb24odHlwZSwgbGlzdGVuZXIpIHtcbiAgaWYgKCFpc0Z1bmN0aW9uKGxpc3RlbmVyKSlcbiAgICB0aHJvdyBUeXBlRXJyb3IoJ2xpc3RlbmVyIG11c3QgYmUgYSBmdW5jdGlvbicpO1xuXG4gIHZhciBmaXJlZCA9IGZhbHNlO1xuXG4gIGZ1bmN0aW9uIGcoKSB7XG4gICAgdGhpcy5yZW1vdmVMaXN0ZW5lcih0eXBlLCBnKTtcblxuICAgIGlmICghZmlyZWQpIHtcbiAgICAgIGZpcmVkID0gdHJ1ZTtcbiAgICAgIGxpc3RlbmVyLmFwcGx5KHRoaXMsIGFyZ3VtZW50cyk7XG4gICAgfVxuICB9XG5cbiAgZy5saXN0ZW5lciA9IGxpc3RlbmVyO1xuICB0aGlzLm9uKHR5cGUsIGcpO1xuXG4gIHJldHVybiB0aGlzO1xufTtcblxuLy8gZW1pdHMgYSAncmVtb3ZlTGlzdGVuZXInIGV2ZW50IGlmZiB0aGUgbGlzdGVuZXIgd2FzIHJlbW92ZWRcbkV2ZW50RW1pdHRlci5wcm90b3R5cGUucmVtb3ZlTGlzdGVuZXIgPSBmdW5jdGlvbih0eXBlLCBsaXN0ZW5lcikge1xuICB2YXIgbGlzdCwgcG9zaXRpb24sIGxlbmd0aCwgaTtcblxuICBpZiAoIWlzRnVuY3Rpb24obGlzdGVuZXIpKVxuICAgIHRocm93IFR5cGVFcnJvcignbGlzdGVuZXIgbXVzdCBiZSBhIGZ1bmN0aW9uJyk7XG5cbiAgaWYgKCF0aGlzLl9ldmVudHMgfHwgIXRoaXMuX2V2ZW50c1t0eXBlXSlcbiAgICByZXR1cm4gdGhpcztcblxuICBsaXN0ID0gdGhpcy5fZXZlbnRzW3R5cGVdO1xuICBsZW5ndGggPSBsaXN0Lmxlbmd0aDtcbiAgcG9zaXRpb24gPSAtMTtcblxuICBpZiAobGlzdCA9PT0gbGlzdGVuZXIgfHxcbiAgICAgIChpc0Z1bmN0aW9uKGxpc3QubGlzdGVuZXIpICYmIGxpc3QubGlzdGVuZXIgPT09IGxpc3RlbmVyKSkge1xuICAgIGRlbGV0ZSB0aGlzLl9ldmVudHNbdHlwZV07XG4gICAgaWYgKHRoaXMuX2V2ZW50cy5yZW1vdmVMaXN0ZW5lcilcbiAgICAgIHRoaXMuZW1pdCgncmVtb3ZlTGlzdGVuZXInLCB0eXBlLCBsaXN0ZW5lcik7XG5cbiAgfSBlbHNlIGlmIChpc09iamVjdChsaXN0KSkge1xuICAgIGZvciAoaSA9IGxlbmd0aDsgaS0tID4gMDspIHtcbiAgICAgIGlmIChsaXN0W2ldID09PSBsaXN0ZW5lciB8fFxuICAgICAgICAgIChsaXN0W2ldLmxpc3RlbmVyICYmIGxpc3RbaV0ubGlzdGVuZXIgPT09IGxpc3RlbmVyKSkge1xuICAgICAgICBwb3NpdGlvbiA9IGk7XG4gICAgICAgIGJyZWFrO1xuICAgICAgfVxuICAgIH1cblxuICAgIGlmIChwb3NpdGlvbiA8IDApXG4gICAgICByZXR1cm4gdGhpcztcblxuICAgIGlmIChsaXN0Lmxlbmd0aCA9PT0gMSkge1xuICAgICAgbGlzdC5sZW5ndGggPSAwO1xuICAgICAgZGVsZXRlIHRoaXMuX2V2ZW50c1t0eXBlXTtcbiAgICB9IGVsc2Uge1xuICAgICAgbGlzdC5zcGxpY2UocG9zaXRpb24sIDEpO1xuICAgIH1cblxuICAgIGlmICh0aGlzLl9ldmVudHMucmVtb3ZlTGlzdGVuZXIpXG4gICAgICB0aGlzLmVtaXQoJ3JlbW92ZUxpc3RlbmVyJywgdHlwZSwgbGlzdGVuZXIpO1xuICB9XG5cbiAgcmV0dXJuIHRoaXM7XG59O1xuXG5FdmVudEVtaXR0ZXIucHJvdG90eXBlLnJlbW92ZUFsbExpc3RlbmVycyA9IGZ1bmN0aW9uKHR5cGUpIHtcbiAgdmFyIGtleSwgbGlzdGVuZXJzO1xuXG4gIGlmICghdGhpcy5fZXZlbnRzKVxuICAgIHJldHVybiB0aGlzO1xuXG4gIC8vIG5vdCBsaXN0ZW5pbmcgZm9yIHJlbW92ZUxpc3RlbmVyLCBubyBuZWVkIHRvIGVtaXRcbiAgaWYgKCF0aGlzLl9ldmVudHMucmVtb3ZlTGlzdGVuZXIpIHtcbiAgICBpZiAoYXJndW1lbnRzLmxlbmd0aCA9PT0gMClcbiAgICAgIHRoaXMuX2V2ZW50cyA9IHt9O1xuICAgIGVsc2UgaWYgKHRoaXMuX2V2ZW50c1t0eXBlXSlcbiAgICAgIGRlbGV0ZSB0aGlzLl9ldmVudHNbdHlwZV07XG4gICAgcmV0dXJuIHRoaXM7XG4gIH1cblxuICAvLyBlbWl0IHJlbW92ZUxpc3RlbmVyIGZvciBhbGwgbGlzdGVuZXJzIG9uIGFsbCBldmVudHNcbiAgaWYgKGFyZ3VtZW50cy5sZW5ndGggPT09IDApIHtcbiAgICBmb3IgKGtleSBpbiB0aGlzLl9ldmVudHMpIHtcbiAgICAgIGlmIChrZXkgPT09ICdyZW1vdmVMaXN0ZW5lcicpIGNvbnRpbnVlO1xuICAgICAgdGhpcy5yZW1vdmVBbGxMaXN0ZW5lcnMoa2V5KTtcbiAgICB9XG4gICAgdGhpcy5yZW1vdmVBbGxMaXN0ZW5lcnMoJ3JlbW92ZUxpc3RlbmVyJyk7XG4gICAgdGhpcy5fZXZlbnRzID0ge307XG4gICAgcmV0dXJuIHRoaXM7XG4gIH1cblxuICBsaXN0ZW5lcnMgPSB0aGlzLl9ldmVudHNbdHlwZV07XG5cbiAgaWYgKGlzRnVuY3Rpb24obGlzdGVuZXJzKSkge1xuICAgIHRoaXMucmVtb3ZlTGlzdGVuZXIodHlwZSwgbGlzdGVuZXJzKTtcbiAgfSBlbHNlIHtcbiAgICAvLyBMSUZPIG9yZGVyXG4gICAgd2hpbGUgKGxpc3RlbmVycy5sZW5ndGgpXG4gICAgICB0aGlzLnJlbW92ZUxpc3RlbmVyKHR5cGUsIGxpc3RlbmVyc1tsaXN0ZW5lcnMubGVuZ3RoIC0gMV0pO1xuICB9XG4gIGRlbGV0ZSB0aGlzLl9ldmVudHNbdHlwZV07XG5cbiAgcmV0dXJuIHRoaXM7XG59O1xuXG5FdmVudEVtaXR0ZXIucHJvdG90eXBlLmxpc3RlbmVycyA9IGZ1bmN0aW9uKHR5cGUpIHtcbiAgdmFyIHJldDtcbiAgaWYgKCF0aGlzLl9ldmVudHMgfHwgIXRoaXMuX2V2ZW50c1t0eXBlXSlcbiAgICByZXQgPSBbXTtcbiAgZWxzZSBpZiAoaXNGdW5jdGlvbih0aGlzLl9ldmVudHNbdHlwZV0pKVxuICAgIHJldCA9IFt0aGlzLl9ldmVudHNbdHlwZV1dO1xuICBlbHNlXG4gICAgcmV0ID0gdGhpcy5fZXZlbnRzW3R5cGVdLnNsaWNlKCk7XG4gIHJldHVybiByZXQ7XG59O1xuXG5FdmVudEVtaXR0ZXIubGlzdGVuZXJDb3VudCA9IGZ1bmN0aW9uKGVtaXR0ZXIsIHR5cGUpIHtcbiAgdmFyIHJldDtcbiAgaWYgKCFlbWl0dGVyLl9ldmVudHMgfHwgIWVtaXR0ZXIuX2V2ZW50c1t0eXBlXSlcbiAgICByZXQgPSAwO1xuICBlbHNlIGlmIChpc0Z1bmN0aW9uKGVtaXR0ZXIuX2V2ZW50c1t0eXBlXSkpXG4gICAgcmV0ID0gMTtcbiAgZWxzZVxuICAgIHJldCA9IGVtaXR0ZXIuX2V2ZW50c1t0eXBlXS5sZW5ndGg7XG4gIHJldHVybiByZXQ7XG59O1xuXG5mdW5jdGlvbiBpc0Z1bmN0aW9uKGFyZykge1xuICByZXR1cm4gdHlwZW9mIGFyZyA9PT0gJ2Z1bmN0aW9uJztcbn1cblxuZnVuY3Rpb24gaXNOdW1iZXIoYXJnKSB7XG4gIHJldHVybiB0eXBlb2YgYXJnID09PSAnbnVtYmVyJztcbn1cblxuZnVuY3Rpb24gaXNPYmplY3QoYXJnKSB7XG4gIHJldHVybiB0eXBlb2YgYXJnID09PSAnb2JqZWN0JyAmJiBhcmcgIT09IG51bGw7XG59XG5cbmZ1bmN0aW9uIGlzVW5kZWZpbmVkKGFyZykge1xuICByZXR1cm4gYXJnID09PSB2b2lkIDA7XG59XG4iLCJpZiAodHlwZW9mIE9iamVjdC5jcmVhdGUgPT09ICdmdW5jdGlvbicpIHtcbiAgLy8gaW1wbGVtZW50YXRpb24gZnJvbSBzdGFuZGFyZCBub2RlLmpzICd1dGlsJyBtb2R1bGVcbiAgbW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbiBpbmhlcml0cyhjdG9yLCBzdXBlckN0b3IpIHtcbiAgICBjdG9yLnN1cGVyXyA9IHN1cGVyQ3RvclxuICAgIGN0b3IucHJvdG90eXBlID0gT2JqZWN0LmNyZWF0ZShzdXBlckN0b3IucHJvdG90eXBlLCB7XG4gICAgICBjb25zdHJ1Y3Rvcjoge1xuICAgICAgICB2YWx1ZTogY3RvcixcbiAgICAgICAgZW51bWVyYWJsZTogZmFsc2UsXG4gICAgICAgIHdyaXRhYmxlOiB0cnVlLFxuICAgICAgICBjb25maWd1cmFibGU6IHRydWVcbiAgICAgIH1cbiAgICB9KTtcbiAgfTtcbn0gZWxzZSB7XG4gIC8vIG9sZCBzY2hvb2wgc2hpbSBmb3Igb2xkIGJyb3dzZXJzXG4gIG1vZHVsZS5leHBvcnRzID0gZnVuY3Rpb24gaW5oZXJpdHMoY3Rvciwgc3VwZXJDdG9yKSB7XG4gICAgY3Rvci5zdXBlcl8gPSBzdXBlckN0b3JcbiAgICB2YXIgVGVtcEN0b3IgPSBmdW5jdGlvbiAoKSB7fVxuICAgIFRlbXBDdG9yLnByb3RvdHlwZSA9IHN1cGVyQ3Rvci5wcm90b3R5cGVcbiAgICBjdG9yLnByb3RvdHlwZSA9IG5ldyBUZW1wQ3RvcigpXG4gICAgY3Rvci5wcm90b3R5cGUuY29uc3RydWN0b3IgPSBjdG9yXG4gIH1cbn1cbiIsIi8vIHNoaW0gZm9yIHVzaW5nIHByb2Nlc3MgaW4gYnJvd3NlclxuXG52YXIgcHJvY2VzcyA9IG1vZHVsZS5leHBvcnRzID0ge307XG5cbnByb2Nlc3MubmV4dFRpY2sgPSAoZnVuY3Rpb24gKCkge1xuICAgIHZhciBjYW5TZXRJbW1lZGlhdGUgPSB0eXBlb2Ygd2luZG93ICE9PSAndW5kZWZpbmVkJ1xuICAgICYmIHdpbmRvdy5zZXRJbW1lZGlhdGU7XG4gICAgdmFyIGNhblBvc3QgPSB0eXBlb2Ygd2luZG93ICE9PSAndW5kZWZpbmVkJ1xuICAgICYmIHdpbmRvdy5wb3N0TWVzc2FnZSAmJiB3aW5kb3cuYWRkRXZlbnRMaXN0ZW5lclxuICAgIDtcblxuICAgIGlmIChjYW5TZXRJbW1lZGlhdGUpIHtcbiAgICAgICAgcmV0dXJuIGZ1bmN0aW9uIChmKSB7IHJldHVybiB3aW5kb3cuc2V0SW1tZWRpYXRlKGYpIH07XG4gICAgfVxuXG4gICAgaWYgKGNhblBvc3QpIHtcbiAgICAgICAgdmFyIHF1ZXVlID0gW107XG4gICAgICAgIHdpbmRvdy5hZGRFdmVudExpc3RlbmVyKCdtZXNzYWdlJywgZnVuY3Rpb24gKGV2KSB7XG4gICAgICAgICAgICB2YXIgc291cmNlID0gZXYuc291cmNlO1xuICAgICAgICAgICAgaWYgKChzb3VyY2UgPT09IHdpbmRvdyB8fCBzb3VyY2UgPT09IG51bGwpICYmIGV2LmRhdGEgPT09ICdwcm9jZXNzLXRpY2snKSB7XG4gICAgICAgICAgICAgICAgZXYuc3RvcFByb3BhZ2F0aW9uKCk7XG4gICAgICAgICAgICAgICAgaWYgKHF1ZXVlLmxlbmd0aCA+IDApIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGZuID0gcXVldWUuc2hpZnQoKTtcbiAgICAgICAgICAgICAgICAgICAgZm4oKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH0sIHRydWUpO1xuXG4gICAgICAgIHJldHVybiBmdW5jdGlvbiBuZXh0VGljayhmbikge1xuICAgICAgICAgICAgcXVldWUucHVzaChmbik7XG4gICAgICAgICAgICB3aW5kb3cucG9zdE1lc3NhZ2UoJ3Byb2Nlc3MtdGljaycsICcqJyk7XG4gICAgICAgIH07XG4gICAgfVxuXG4gICAgcmV0dXJuIGZ1bmN0aW9uIG5leHRUaWNrKGZuKSB7XG4gICAgICAgIHNldFRpbWVvdXQoZm4sIDApO1xuICAgIH07XG59KSgpO1xuXG5wcm9jZXNzLnRpdGxlID0gJ2Jyb3dzZXInO1xucHJvY2Vzcy5icm93c2VyID0gdHJ1ZTtcbnByb2Nlc3MuZW52ID0ge307XG5wcm9jZXNzLmFyZ3YgPSBbXTtcblxuZnVuY3Rpb24gbm9vcCgpIHt9XG5cbnByb2Nlc3Mub24gPSBub29wO1xucHJvY2Vzcy5hZGRMaXN0ZW5lciA9IG5vb3A7XG5wcm9jZXNzLm9uY2UgPSBub29wO1xucHJvY2Vzcy5vZmYgPSBub29wO1xucHJvY2Vzcy5yZW1vdmVMaXN0ZW5lciA9IG5vb3A7XG5wcm9jZXNzLnJlbW92ZUFsbExpc3RlbmVycyA9IG5vb3A7XG5wcm9jZXNzLmVtaXQgPSBub29wO1xuXG5wcm9jZXNzLmJpbmRpbmcgPSBmdW5jdGlvbiAobmFtZSkge1xuICAgIHRocm93IG5ldyBFcnJvcigncHJvY2Vzcy5iaW5kaW5nIGlzIG5vdCBzdXBwb3J0ZWQnKTtcbn1cblxuLy8gVE9ETyhzaHR5bG1hbilcbnByb2Nlc3MuY3dkID0gZnVuY3Rpb24gKCkgeyByZXR1cm4gJy8nIH07XG5wcm9jZXNzLmNoZGlyID0gZnVuY3Rpb24gKGRpcikge1xuICAgIHRocm93IG5ldyBFcnJvcigncHJvY2Vzcy5jaGRpciBpcyBub3Qgc3VwcG9ydGVkJyk7XG59O1xuIiwibW9kdWxlLmV4cG9ydHMgPSBmdW5jdGlvbiBpc0J1ZmZlcihhcmcpIHtcbiAgcmV0dXJuIGFyZyAmJiB0eXBlb2YgYXJnID09PSAnb2JqZWN0J1xuICAgICYmIHR5cGVvZiBhcmcuY29weSA9PT0gJ2Z1bmN0aW9uJ1xuICAgICYmIHR5cGVvZiBhcmcuZmlsbCA9PT0gJ2Z1bmN0aW9uJ1xuICAgICYmIHR5cGVvZiBhcmcucmVhZFVJbnQ4ID09PSAnZnVuY3Rpb24nO1xufSIsIihmdW5jdGlvbiAocHJvY2VzcyxnbG9iYWwpe1xuLy8gQ29weXJpZ2h0IEpveWVudCwgSW5jLiBhbmQgb3RoZXIgTm9kZSBjb250cmlidXRvcnMuXG4vL1xuLy8gUGVybWlzc2lvbiBpcyBoZXJlYnkgZ3JhbnRlZCwgZnJlZSBvZiBjaGFyZ2UsIHRvIGFueSBwZXJzb24gb2J0YWluaW5nIGFcbi8vIGNvcHkgb2YgdGhpcyBzb2Z0d2FyZSBhbmQgYXNzb2NpYXRlZCBkb2N1bWVudGF0aW9uIGZpbGVzICh0aGVcbi8vIFwiU29mdHdhcmVcIiksIHRvIGRlYWwgaW4gdGhlIFNvZnR3YXJlIHdpdGhvdXQgcmVzdHJpY3Rpb24sIGluY2x1ZGluZ1xuLy8gd2l0aG91dCBsaW1pdGF0aW9uIHRoZSByaWdodHMgdG8gdXNlLCBjb3B5LCBtb2RpZnksIG1lcmdlLCBwdWJsaXNoLFxuLy8gZGlzdHJpYnV0ZSwgc3VibGljZW5zZSwgYW5kL29yIHNlbGwgY29waWVzIG9mIHRoZSBTb2Z0d2FyZSwgYW5kIHRvIHBlcm1pdFxuLy8gcGVyc29ucyB0byB3aG9tIHRoZSBTb2Z0d2FyZSBpcyBmdXJuaXNoZWQgdG8gZG8gc28sIHN1YmplY3QgdG8gdGhlXG4vLyBmb2xsb3dpbmcgY29uZGl0aW9uczpcbi8vXG4vLyBUaGUgYWJvdmUgY29weXJpZ2h0IG5vdGljZSBhbmQgdGhpcyBwZXJtaXNzaW9uIG5vdGljZSBzaGFsbCBiZSBpbmNsdWRlZFxuLy8gaW4gYWxsIGNvcGllcyBvciBzdWJzdGFudGlhbCBwb3J0aW9ucyBvZiB0aGUgU29mdHdhcmUuXG4vL1xuLy8gVEhFIFNPRlRXQVJFIElTIFBST1ZJREVEIFwiQVMgSVNcIiwgV0lUSE9VVCBXQVJSQU5UWSBPRiBBTlkgS0lORCwgRVhQUkVTU1xuLy8gT1IgSU1QTElFRCwgSU5DTFVESU5HIEJVVCBOT1QgTElNSVRFRCBUTyBUSEUgV0FSUkFOVElFUyBPRlxuLy8gTUVSQ0hBTlRBQklMSVRZLCBGSVRORVNTIEZPUiBBIFBBUlRJQ1VMQVIgUFVSUE9TRSBBTkQgTk9OSU5GUklOR0VNRU5ULiBJTlxuLy8gTk8gRVZFTlQgU0hBTEwgVEhFIEFVVEhPUlMgT1IgQ09QWVJJR0hUIEhPTERFUlMgQkUgTElBQkxFIEZPUiBBTlkgQ0xBSU0sXG4vLyBEQU1BR0VTIE9SIE9USEVSIExJQUJJTElUWSwgV0hFVEhFUiBJTiBBTiBBQ1RJT04gT0YgQ09OVFJBQ1QsIFRPUlQgT1Jcbi8vIE9USEVSV0lTRSwgQVJJU0lORyBGUk9NLCBPVVQgT0YgT1IgSU4gQ09OTkVDVElPTiBXSVRIIFRIRSBTT0ZUV0FSRSBPUiBUSEVcbi8vIFVTRSBPUiBPVEhFUiBERUFMSU5HUyBJTiBUSEUgU09GVFdBUkUuXG5cbnZhciBmb3JtYXRSZWdFeHAgPSAvJVtzZGolXS9nO1xuZXhwb3J0cy5mb3JtYXQgPSBmdW5jdGlvbihmKSB7XG4gIGlmICghaXNTdHJpbmcoZikpIHtcbiAgICB2YXIgb2JqZWN0cyA9IFtdO1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgYXJndW1lbnRzLmxlbmd0aDsgaSsrKSB7XG4gICAgICBvYmplY3RzLnB1c2goaW5zcGVjdChhcmd1bWVudHNbaV0pKTtcbiAgICB9XG4gICAgcmV0dXJuIG9iamVjdHMuam9pbignICcpO1xuICB9XG5cbiAgdmFyIGkgPSAxO1xuICB2YXIgYXJncyA9IGFyZ3VtZW50cztcbiAgdmFyIGxlbiA9IGFyZ3MubGVuZ3RoO1xuICB2YXIgc3RyID0gU3RyaW5nKGYpLnJlcGxhY2UoZm9ybWF0UmVnRXhwLCBmdW5jdGlvbih4KSB7XG4gICAgaWYgKHggPT09ICclJScpIHJldHVybiAnJSc7XG4gICAgaWYgKGkgPj0gbGVuKSByZXR1cm4geDtcbiAgICBzd2l0Y2ggKHgpIHtcbiAgICAgIGNhc2UgJyVzJzogcmV0dXJuIFN0cmluZyhhcmdzW2krK10pO1xuICAgICAgY2FzZSAnJWQnOiByZXR1cm4gTnVtYmVyKGFyZ3NbaSsrXSk7XG4gICAgICBjYXNlICclaic6XG4gICAgICAgIHRyeSB7XG4gICAgICAgICAgcmV0dXJuIEpTT04uc3RyaW5naWZ5KGFyZ3NbaSsrXSk7XG4gICAgICAgIH0gY2F0Y2ggKF8pIHtcbiAgICAgICAgICByZXR1cm4gJ1tDaXJjdWxhcl0nO1xuICAgICAgICB9XG4gICAgICBkZWZhdWx0OlxuICAgICAgICByZXR1cm4geDtcbiAgICB9XG4gIH0pO1xuICBmb3IgKHZhciB4ID0gYXJnc1tpXTsgaSA8IGxlbjsgeCA9IGFyZ3NbKytpXSkge1xuICAgIGlmIChpc051bGwoeCkgfHwgIWlzT2JqZWN0KHgpKSB7XG4gICAgICBzdHIgKz0gJyAnICsgeDtcbiAgICB9IGVsc2Uge1xuICAgICAgc3RyICs9ICcgJyArIGluc3BlY3QoeCk7XG4gICAgfVxuICB9XG4gIHJldHVybiBzdHI7XG59O1xuXG5cbi8vIE1hcmsgdGhhdCBhIG1ldGhvZCBzaG91bGQgbm90IGJlIHVzZWQuXG4vLyBSZXR1cm5zIGEgbW9kaWZpZWQgZnVuY3Rpb24gd2hpY2ggd2FybnMgb25jZSBieSBkZWZhdWx0LlxuLy8gSWYgLS1uby1kZXByZWNhdGlvbiBpcyBzZXQsIHRoZW4gaXQgaXMgYSBuby1vcC5cbmV4cG9ydHMuZGVwcmVjYXRlID0gZnVuY3Rpb24oZm4sIG1zZykge1xuICAvLyBBbGxvdyBmb3IgZGVwcmVjYXRpbmcgdGhpbmdzIGluIHRoZSBwcm9jZXNzIG9mIHN0YXJ0aW5nIHVwLlxuICBpZiAoaXNVbmRlZmluZWQoZ2xvYmFsLnByb2Nlc3MpKSB7XG4gICAgcmV0dXJuIGZ1bmN0aW9uKCkge1xuICAgICAgcmV0dXJuIGV4cG9ydHMuZGVwcmVjYXRlKGZuLCBtc2cpLmFwcGx5KHRoaXMsIGFyZ3VtZW50cyk7XG4gICAgfTtcbiAgfVxuXG4gIGlmIChwcm9jZXNzLm5vRGVwcmVjYXRpb24gPT09IHRydWUpIHtcbiAgICByZXR1cm4gZm47XG4gIH1cblxuICB2YXIgd2FybmVkID0gZmFsc2U7XG4gIGZ1bmN0aW9uIGRlcHJlY2F0ZWQoKSB7XG4gICAgaWYgKCF3YXJuZWQpIHtcbiAgICAgIGlmIChwcm9jZXNzLnRocm93RGVwcmVjYXRpb24pIHtcbiAgICAgICAgdGhyb3cgbmV3IEVycm9yKG1zZyk7XG4gICAgICB9IGVsc2UgaWYgKHByb2Nlc3MudHJhY2VEZXByZWNhdGlvbikge1xuICAgICAgICBjb25zb2xlLnRyYWNlKG1zZyk7XG4gICAgICB9IGVsc2Uge1xuICAgICAgICBjb25zb2xlLmVycm9yKG1zZyk7XG4gICAgICB9XG4gICAgICB3YXJuZWQgPSB0cnVlO1xuICAgIH1cbiAgICByZXR1cm4gZm4uYXBwbHkodGhpcywgYXJndW1lbnRzKTtcbiAgfVxuXG4gIHJldHVybiBkZXByZWNhdGVkO1xufTtcblxuXG52YXIgZGVidWdzID0ge307XG52YXIgZGVidWdFbnZpcm9uO1xuZXhwb3J0cy5kZWJ1Z2xvZyA9IGZ1bmN0aW9uKHNldCkge1xuICBpZiAoaXNVbmRlZmluZWQoZGVidWdFbnZpcm9uKSlcbiAgICBkZWJ1Z0Vudmlyb24gPSBwcm9jZXNzLmVudi5OT0RFX0RFQlVHIHx8ICcnO1xuICBzZXQgPSBzZXQudG9VcHBlckNhc2UoKTtcbiAgaWYgKCFkZWJ1Z3Nbc2V0XSkge1xuICAgIGlmIChuZXcgUmVnRXhwKCdcXFxcYicgKyBzZXQgKyAnXFxcXGInLCAnaScpLnRlc3QoZGVidWdFbnZpcm9uKSkge1xuICAgICAgdmFyIHBpZCA9IHByb2Nlc3MucGlkO1xuICAgICAgZGVidWdzW3NldF0gPSBmdW5jdGlvbigpIHtcbiAgICAgICAgdmFyIG1zZyA9IGV4cG9ydHMuZm9ybWF0LmFwcGx5KGV4cG9ydHMsIGFyZ3VtZW50cyk7XG4gICAgICAgIGNvbnNvbGUuZXJyb3IoJyVzICVkOiAlcycsIHNldCwgcGlkLCBtc2cpO1xuICAgICAgfTtcbiAgICB9IGVsc2Uge1xuICAgICAgZGVidWdzW3NldF0gPSBmdW5jdGlvbigpIHt9O1xuICAgIH1cbiAgfVxuICByZXR1cm4gZGVidWdzW3NldF07XG59O1xuXG5cbi8qKlxuICogRWNob3MgdGhlIHZhbHVlIG9mIGEgdmFsdWUuIFRyeXMgdG8gcHJpbnQgdGhlIHZhbHVlIG91dFxuICogaW4gdGhlIGJlc3Qgd2F5IHBvc3NpYmxlIGdpdmVuIHRoZSBkaWZmZXJlbnQgdHlwZXMuXG4gKlxuICogQHBhcmFtIHtPYmplY3R9IG9iaiBUaGUgb2JqZWN0IHRvIHByaW50IG91dC5cbiAqIEBwYXJhbSB7T2JqZWN0fSBvcHRzIE9wdGlvbmFsIG9wdGlvbnMgb2JqZWN0IHRoYXQgYWx0ZXJzIHRoZSBvdXRwdXQuXG4gKi9cbi8qIGxlZ2FjeTogb2JqLCBzaG93SGlkZGVuLCBkZXB0aCwgY29sb3JzKi9cbmZ1bmN0aW9uIGluc3BlY3Qob2JqLCBvcHRzKSB7XG4gIC8vIGRlZmF1bHQgb3B0aW9uc1xuICB2YXIgY3R4ID0ge1xuICAgIHNlZW46IFtdLFxuICAgIHN0eWxpemU6IHN0eWxpemVOb0NvbG9yXG4gIH07XG4gIC8vIGxlZ2FjeS4uLlxuICBpZiAoYXJndW1lbnRzLmxlbmd0aCA+PSAzKSBjdHguZGVwdGggPSBhcmd1bWVudHNbMl07XG4gIGlmIChhcmd1bWVudHMubGVuZ3RoID49IDQpIGN0eC5jb2xvcnMgPSBhcmd1bWVudHNbM107XG4gIGlmIChpc0Jvb2xlYW4ob3B0cykpIHtcbiAgICAvLyBsZWdhY3kuLi5cbiAgICBjdHguc2hvd0hpZGRlbiA9IG9wdHM7XG4gIH0gZWxzZSBpZiAob3B0cykge1xuICAgIC8vIGdvdCBhbiBcIm9wdGlvbnNcIiBvYmplY3RcbiAgICBleHBvcnRzLl9leHRlbmQoY3R4LCBvcHRzKTtcbiAgfVxuICAvLyBzZXQgZGVmYXVsdCBvcHRpb25zXG4gIGlmIChpc1VuZGVmaW5lZChjdHguc2hvd0hpZGRlbikpIGN0eC5zaG93SGlkZGVuID0gZmFsc2U7XG4gIGlmIChpc1VuZGVmaW5lZChjdHguZGVwdGgpKSBjdHguZGVwdGggPSAyO1xuICBpZiAoaXNVbmRlZmluZWQoY3R4LmNvbG9ycykpIGN0eC5jb2xvcnMgPSBmYWxzZTtcbiAgaWYgKGlzVW5kZWZpbmVkKGN0eC5jdXN0b21JbnNwZWN0KSkgY3R4LmN1c3RvbUluc3BlY3QgPSB0cnVlO1xuICBpZiAoY3R4LmNvbG9ycykgY3R4LnN0eWxpemUgPSBzdHlsaXplV2l0aENvbG9yO1xuICByZXR1cm4gZm9ybWF0VmFsdWUoY3R4LCBvYmosIGN0eC5kZXB0aCk7XG59XG5leHBvcnRzLmluc3BlY3QgPSBpbnNwZWN0O1xuXG5cbi8vIGh0dHA6Ly9lbi53aWtpcGVkaWEub3JnL3dpa2kvQU5TSV9lc2NhcGVfY29kZSNncmFwaGljc1xuaW5zcGVjdC5jb2xvcnMgPSB7XG4gICdib2xkJyA6IFsxLCAyMl0sXG4gICdpdGFsaWMnIDogWzMsIDIzXSxcbiAgJ3VuZGVybGluZScgOiBbNCwgMjRdLFxuICAnaW52ZXJzZScgOiBbNywgMjddLFxuICAnd2hpdGUnIDogWzM3LCAzOV0sXG4gICdncmV5JyA6IFs5MCwgMzldLFxuICAnYmxhY2snIDogWzMwLCAzOV0sXG4gICdibHVlJyA6IFszNCwgMzldLFxuICAnY3lhbicgOiBbMzYsIDM5XSxcbiAgJ2dyZWVuJyA6IFszMiwgMzldLFxuICAnbWFnZW50YScgOiBbMzUsIDM5XSxcbiAgJ3JlZCcgOiBbMzEsIDM5XSxcbiAgJ3llbGxvdycgOiBbMzMsIDM5XVxufTtcblxuLy8gRG9uJ3QgdXNlICdibHVlJyBub3QgdmlzaWJsZSBvbiBjbWQuZXhlXG5pbnNwZWN0LnN0eWxlcyA9IHtcbiAgJ3NwZWNpYWwnOiAnY3lhbicsXG4gICdudW1iZXInOiAneWVsbG93JyxcbiAgJ2Jvb2xlYW4nOiAneWVsbG93JyxcbiAgJ3VuZGVmaW5lZCc6ICdncmV5JyxcbiAgJ251bGwnOiAnYm9sZCcsXG4gICdzdHJpbmcnOiAnZ3JlZW4nLFxuICAnZGF0ZSc6ICdtYWdlbnRhJyxcbiAgLy8gXCJuYW1lXCI6IGludGVudGlvbmFsbHkgbm90IHN0eWxpbmdcbiAgJ3JlZ2V4cCc6ICdyZWQnXG59O1xuXG5cbmZ1bmN0aW9uIHN0eWxpemVXaXRoQ29sb3Ioc3RyLCBzdHlsZVR5cGUpIHtcbiAgdmFyIHN0eWxlID0gaW5zcGVjdC5zdHlsZXNbc3R5bGVUeXBlXTtcblxuICBpZiAoc3R5bGUpIHtcbiAgICByZXR1cm4gJ1xcdTAwMWJbJyArIGluc3BlY3QuY29sb3JzW3N0eWxlXVswXSArICdtJyArIHN0ciArXG4gICAgICAgICAgICdcXHUwMDFiWycgKyBpbnNwZWN0LmNvbG9yc1tzdHlsZV1bMV0gKyAnbSc7XG4gIH0gZWxzZSB7XG4gICAgcmV0dXJuIHN0cjtcbiAgfVxufVxuXG5cbmZ1bmN0aW9uIHN0eWxpemVOb0NvbG9yKHN0ciwgc3R5bGVUeXBlKSB7XG4gIHJldHVybiBzdHI7XG59XG5cblxuZnVuY3Rpb24gYXJyYXlUb0hhc2goYXJyYXkpIHtcbiAgdmFyIGhhc2ggPSB7fTtcblxuICBhcnJheS5mb3JFYWNoKGZ1bmN0aW9uKHZhbCwgaWR4KSB7XG4gICAgaGFzaFt2YWxdID0gdHJ1ZTtcbiAgfSk7XG5cbiAgcmV0dXJuIGhhc2g7XG59XG5cblxuZnVuY3Rpb24gZm9ybWF0VmFsdWUoY3R4LCB2YWx1ZSwgcmVjdXJzZVRpbWVzKSB7XG4gIC8vIFByb3ZpZGUgYSBob29rIGZvciB1c2VyLXNwZWNpZmllZCBpbnNwZWN0IGZ1bmN0aW9ucy5cbiAgLy8gQ2hlY2sgdGhhdCB2YWx1ZSBpcyBhbiBvYmplY3Qgd2l0aCBhbiBpbnNwZWN0IGZ1bmN0aW9uIG9uIGl0XG4gIGlmIChjdHguY3VzdG9tSW5zcGVjdCAmJlxuICAgICAgdmFsdWUgJiZcbiAgICAgIGlzRnVuY3Rpb24odmFsdWUuaW5zcGVjdCkgJiZcbiAgICAgIC8vIEZpbHRlciBvdXQgdGhlIHV0aWwgbW9kdWxlLCBpdCdzIGluc3BlY3QgZnVuY3Rpb24gaXMgc3BlY2lhbFxuICAgICAgdmFsdWUuaW5zcGVjdCAhPT0gZXhwb3J0cy5pbnNwZWN0ICYmXG4gICAgICAvLyBBbHNvIGZpbHRlciBvdXQgYW55IHByb3RvdHlwZSBvYmplY3RzIHVzaW5nIHRoZSBjaXJjdWxhciBjaGVjay5cbiAgICAgICEodmFsdWUuY29uc3RydWN0b3IgJiYgdmFsdWUuY29uc3RydWN0b3IucHJvdG90eXBlID09PSB2YWx1ZSkpIHtcbiAgICB2YXIgcmV0ID0gdmFsdWUuaW5zcGVjdChyZWN1cnNlVGltZXMsIGN0eCk7XG4gICAgaWYgKCFpc1N0cmluZyhyZXQpKSB7XG4gICAgICByZXQgPSBmb3JtYXRWYWx1ZShjdHgsIHJldCwgcmVjdXJzZVRpbWVzKTtcbiAgICB9XG4gICAgcmV0dXJuIHJldDtcbiAgfVxuXG4gIC8vIFByaW1pdGl2ZSB0eXBlcyBjYW5ub3QgaGF2ZSBwcm9wZXJ0aWVzXG4gIHZhciBwcmltaXRpdmUgPSBmb3JtYXRQcmltaXRpdmUoY3R4LCB2YWx1ZSk7XG4gIGlmIChwcmltaXRpdmUpIHtcbiAgICByZXR1cm4gcHJpbWl0aXZlO1xuICB9XG5cbiAgLy8gTG9vayB1cCB0aGUga2V5cyBvZiB0aGUgb2JqZWN0LlxuICB2YXIga2V5cyA9IE9iamVjdC5rZXlzKHZhbHVlKTtcbiAgdmFyIHZpc2libGVLZXlzID0gYXJyYXlUb0hhc2goa2V5cyk7XG5cbiAgaWYgKGN0eC5zaG93SGlkZGVuKSB7XG4gICAga2V5cyA9IE9iamVjdC5nZXRPd25Qcm9wZXJ0eU5hbWVzKHZhbHVlKTtcbiAgfVxuXG4gIC8vIElFIGRvZXNuJ3QgbWFrZSBlcnJvciBmaWVsZHMgbm9uLWVudW1lcmFibGVcbiAgLy8gaHR0cDovL21zZG4ubWljcm9zb2Z0LmNvbS9lbi11cy9saWJyYXJ5L2llL2R3dzUyc2J0KHY9dnMuOTQpLmFzcHhcbiAgaWYgKGlzRXJyb3IodmFsdWUpXG4gICAgICAmJiAoa2V5cy5pbmRleE9mKCdtZXNzYWdlJykgPj0gMCB8fCBrZXlzLmluZGV4T2YoJ2Rlc2NyaXB0aW9uJykgPj0gMCkpIHtcbiAgICByZXR1cm4gZm9ybWF0RXJyb3IodmFsdWUpO1xuICB9XG5cbiAgLy8gU29tZSB0eXBlIG9mIG9iamVjdCB3aXRob3V0IHByb3BlcnRpZXMgY2FuIGJlIHNob3J0Y3V0dGVkLlxuICBpZiAoa2V5cy5sZW5ndGggPT09IDApIHtcbiAgICBpZiAoaXNGdW5jdGlvbih2YWx1ZSkpIHtcbiAgICAgIHZhciBuYW1lID0gdmFsdWUubmFtZSA/ICc6ICcgKyB2YWx1ZS5uYW1lIDogJyc7XG4gICAgICByZXR1cm4gY3R4LnN0eWxpemUoJ1tGdW5jdGlvbicgKyBuYW1lICsgJ10nLCAnc3BlY2lhbCcpO1xuICAgIH1cbiAgICBpZiAoaXNSZWdFeHAodmFsdWUpKSB7XG4gICAgICByZXR1cm4gY3R4LnN0eWxpemUoUmVnRXhwLnByb3RvdHlwZS50b1N0cmluZy5jYWxsKHZhbHVlKSwgJ3JlZ2V4cCcpO1xuICAgIH1cbiAgICBpZiAoaXNEYXRlKHZhbHVlKSkge1xuICAgICAgcmV0dXJuIGN0eC5zdHlsaXplKERhdGUucHJvdG90eXBlLnRvU3RyaW5nLmNhbGwodmFsdWUpLCAnZGF0ZScpO1xuICAgIH1cbiAgICBpZiAoaXNFcnJvcih2YWx1ZSkpIHtcbiAgICAgIHJldHVybiBmb3JtYXRFcnJvcih2YWx1ZSk7XG4gICAgfVxuICB9XG5cbiAgdmFyIGJhc2UgPSAnJywgYXJyYXkgPSBmYWxzZSwgYnJhY2VzID0gWyd7JywgJ30nXTtcblxuICAvLyBNYWtlIEFycmF5IHNheSB0aGF0IHRoZXkgYXJlIEFycmF5XG4gIGlmIChpc0FycmF5KHZhbHVlKSkge1xuICAgIGFycmF5ID0gdHJ1ZTtcbiAgICBicmFjZXMgPSBbJ1snLCAnXSddO1xuICB9XG5cbiAgLy8gTWFrZSBmdW5jdGlvbnMgc2F5IHRoYXQgdGhleSBhcmUgZnVuY3Rpb25zXG4gIGlmIChpc0Z1bmN0aW9uKHZhbHVlKSkge1xuICAgIHZhciBuID0gdmFsdWUubmFtZSA/ICc6ICcgKyB2YWx1ZS5uYW1lIDogJyc7XG4gICAgYmFzZSA9ICcgW0Z1bmN0aW9uJyArIG4gKyAnXSc7XG4gIH1cblxuICAvLyBNYWtlIFJlZ0V4cHMgc2F5IHRoYXQgdGhleSBhcmUgUmVnRXhwc1xuICBpZiAoaXNSZWdFeHAodmFsdWUpKSB7XG4gICAgYmFzZSA9ICcgJyArIFJlZ0V4cC5wcm90b3R5cGUudG9TdHJpbmcuY2FsbCh2YWx1ZSk7XG4gIH1cblxuICAvLyBNYWtlIGRhdGVzIHdpdGggcHJvcGVydGllcyBmaXJzdCBzYXkgdGhlIGRhdGVcbiAgaWYgKGlzRGF0ZSh2YWx1ZSkpIHtcbiAgICBiYXNlID0gJyAnICsgRGF0ZS5wcm90b3R5cGUudG9VVENTdHJpbmcuY2FsbCh2YWx1ZSk7XG4gIH1cblxuICAvLyBNYWtlIGVycm9yIHdpdGggbWVzc2FnZSBmaXJzdCBzYXkgdGhlIGVycm9yXG4gIGlmIChpc0Vycm9yKHZhbHVlKSkge1xuICAgIGJhc2UgPSAnICcgKyBmb3JtYXRFcnJvcih2YWx1ZSk7XG4gIH1cblxuICBpZiAoa2V5cy5sZW5ndGggPT09IDAgJiYgKCFhcnJheSB8fCB2YWx1ZS5sZW5ndGggPT0gMCkpIHtcbiAgICByZXR1cm4gYnJhY2VzWzBdICsgYmFzZSArIGJyYWNlc1sxXTtcbiAgfVxuXG4gIGlmIChyZWN1cnNlVGltZXMgPCAwKSB7XG4gICAgaWYgKGlzUmVnRXhwKHZhbHVlKSkge1xuICAgICAgcmV0dXJuIGN0eC5zdHlsaXplKFJlZ0V4cC5wcm90b3R5cGUudG9TdHJpbmcuY2FsbCh2YWx1ZSksICdyZWdleHAnKTtcbiAgICB9IGVsc2Uge1xuICAgICAgcmV0dXJuIGN0eC5zdHlsaXplKCdbT2JqZWN0XScsICdzcGVjaWFsJyk7XG4gICAgfVxuICB9XG5cbiAgY3R4LnNlZW4ucHVzaCh2YWx1ZSk7XG5cbiAgdmFyIG91dHB1dDtcbiAgaWYgKGFycmF5KSB7XG4gICAgb3V0cHV0ID0gZm9ybWF0QXJyYXkoY3R4LCB2YWx1ZSwgcmVjdXJzZVRpbWVzLCB2aXNpYmxlS2V5cywga2V5cyk7XG4gIH0gZWxzZSB7XG4gICAgb3V0cHV0ID0ga2V5cy5tYXAoZnVuY3Rpb24oa2V5KSB7XG4gICAgICByZXR1cm4gZm9ybWF0UHJvcGVydHkoY3R4LCB2YWx1ZSwgcmVjdXJzZVRpbWVzLCB2aXNpYmxlS2V5cywga2V5LCBhcnJheSk7XG4gICAgfSk7XG4gIH1cblxuICBjdHguc2Vlbi5wb3AoKTtcblxuICByZXR1cm4gcmVkdWNlVG9TaW5nbGVTdHJpbmcob3V0cHV0LCBiYXNlLCBicmFjZXMpO1xufVxuXG5cbmZ1bmN0aW9uIGZvcm1hdFByaW1pdGl2ZShjdHgsIHZhbHVlKSB7XG4gIGlmIChpc1VuZGVmaW5lZCh2YWx1ZSkpXG4gICAgcmV0dXJuIGN0eC5zdHlsaXplKCd1bmRlZmluZWQnLCAndW5kZWZpbmVkJyk7XG4gIGlmIChpc1N0cmluZyh2YWx1ZSkpIHtcbiAgICB2YXIgc2ltcGxlID0gJ1xcJycgKyBKU09OLnN0cmluZ2lmeSh2YWx1ZSkucmVwbGFjZSgvXlwifFwiJC9nLCAnJylcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIC5yZXBsYWNlKC8nL2csIFwiXFxcXCdcIilcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIC5yZXBsYWNlKC9cXFxcXCIvZywgJ1wiJykgKyAnXFwnJztcbiAgICByZXR1cm4gY3R4LnN0eWxpemUoc2ltcGxlLCAnc3RyaW5nJyk7XG4gIH1cbiAgaWYgKGlzTnVtYmVyKHZhbHVlKSlcbiAgICByZXR1cm4gY3R4LnN0eWxpemUoJycgKyB2YWx1ZSwgJ251bWJlcicpO1xuICBpZiAoaXNCb29sZWFuKHZhbHVlKSlcbiAgICByZXR1cm4gY3R4LnN0eWxpemUoJycgKyB2YWx1ZSwgJ2Jvb2xlYW4nKTtcbiAgLy8gRm9yIHNvbWUgcmVhc29uIHR5cGVvZiBudWxsIGlzIFwib2JqZWN0XCIsIHNvIHNwZWNpYWwgY2FzZSBoZXJlLlxuICBpZiAoaXNOdWxsKHZhbHVlKSlcbiAgICByZXR1cm4gY3R4LnN0eWxpemUoJ251bGwnLCAnbnVsbCcpO1xufVxuXG5cbmZ1bmN0aW9uIGZvcm1hdEVycm9yKHZhbHVlKSB7XG4gIHJldHVybiAnWycgKyBFcnJvci5wcm90b3R5cGUudG9TdHJpbmcuY2FsbCh2YWx1ZSkgKyAnXSc7XG59XG5cblxuZnVuY3Rpb24gZm9ybWF0QXJyYXkoY3R4LCB2YWx1ZSwgcmVjdXJzZVRpbWVzLCB2aXNpYmxlS2V5cywga2V5cykge1xuICB2YXIgb3V0cHV0ID0gW107XG4gIGZvciAodmFyIGkgPSAwLCBsID0gdmFsdWUubGVuZ3RoOyBpIDwgbDsgKytpKSB7XG4gICAgaWYgKGhhc093blByb3BlcnR5KHZhbHVlLCBTdHJpbmcoaSkpKSB7XG4gICAgICBvdXRwdXQucHVzaChmb3JtYXRQcm9wZXJ0eShjdHgsIHZhbHVlLCByZWN1cnNlVGltZXMsIHZpc2libGVLZXlzLFxuICAgICAgICAgIFN0cmluZyhpKSwgdHJ1ZSkpO1xuICAgIH0gZWxzZSB7XG4gICAgICBvdXRwdXQucHVzaCgnJyk7XG4gICAgfVxuICB9XG4gIGtleXMuZm9yRWFjaChmdW5jdGlvbihrZXkpIHtcbiAgICBpZiAoIWtleS5tYXRjaCgvXlxcZCskLykpIHtcbiAgICAgIG91dHB1dC5wdXNoKGZvcm1hdFByb3BlcnR5KGN0eCwgdmFsdWUsIHJlY3Vyc2VUaW1lcywgdmlzaWJsZUtleXMsXG4gICAgICAgICAga2V5LCB0cnVlKSk7XG4gICAgfVxuICB9KTtcbiAgcmV0dXJuIG91dHB1dDtcbn1cblxuXG5mdW5jdGlvbiBmb3JtYXRQcm9wZXJ0eShjdHgsIHZhbHVlLCByZWN1cnNlVGltZXMsIHZpc2libGVLZXlzLCBrZXksIGFycmF5KSB7XG4gIHZhciBuYW1lLCBzdHIsIGRlc2M7XG4gIGRlc2MgPSBPYmplY3QuZ2V0T3duUHJvcGVydHlEZXNjcmlwdG9yKHZhbHVlLCBrZXkpIHx8IHsgdmFsdWU6IHZhbHVlW2tleV0gfTtcbiAgaWYgKGRlc2MuZ2V0KSB7XG4gICAgaWYgKGRlc2Muc2V0KSB7XG4gICAgICBzdHIgPSBjdHguc3R5bGl6ZSgnW0dldHRlci9TZXR0ZXJdJywgJ3NwZWNpYWwnKTtcbiAgICB9IGVsc2Uge1xuICAgICAgc3RyID0gY3R4LnN0eWxpemUoJ1tHZXR0ZXJdJywgJ3NwZWNpYWwnKTtcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgaWYgKGRlc2Muc2V0KSB7XG4gICAgICBzdHIgPSBjdHguc3R5bGl6ZSgnW1NldHRlcl0nLCAnc3BlY2lhbCcpO1xuICAgIH1cbiAgfVxuICBpZiAoIWhhc093blByb3BlcnR5KHZpc2libGVLZXlzLCBrZXkpKSB7XG4gICAgbmFtZSA9ICdbJyArIGtleSArICddJztcbiAgfVxuICBpZiAoIXN0cikge1xuICAgIGlmIChjdHguc2Vlbi5pbmRleE9mKGRlc2MudmFsdWUpIDwgMCkge1xuICAgICAgaWYgKGlzTnVsbChyZWN1cnNlVGltZXMpKSB7XG4gICAgICAgIHN0ciA9IGZvcm1hdFZhbHVlKGN0eCwgZGVzYy52YWx1ZSwgbnVsbCk7XG4gICAgICB9IGVsc2Uge1xuICAgICAgICBzdHIgPSBmb3JtYXRWYWx1ZShjdHgsIGRlc2MudmFsdWUsIHJlY3Vyc2VUaW1lcyAtIDEpO1xuICAgICAgfVxuICAgICAgaWYgKHN0ci5pbmRleE9mKCdcXG4nKSA+IC0xKSB7XG4gICAgICAgIGlmIChhcnJheSkge1xuICAgICAgICAgIHN0ciA9IHN0ci5zcGxpdCgnXFxuJykubWFwKGZ1bmN0aW9uKGxpbmUpIHtcbiAgICAgICAgICAgIHJldHVybiAnICAnICsgbGluZTtcbiAgICAgICAgICB9KS5qb2luKCdcXG4nKS5zdWJzdHIoMik7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgc3RyID0gJ1xcbicgKyBzdHIuc3BsaXQoJ1xcbicpLm1hcChmdW5jdGlvbihsaW5lKSB7XG4gICAgICAgICAgICByZXR1cm4gJyAgICcgKyBsaW5lO1xuICAgICAgICAgIH0pLmpvaW4oJ1xcbicpO1xuICAgICAgICB9XG4gICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgIHN0ciA9IGN0eC5zdHlsaXplKCdbQ2lyY3VsYXJdJywgJ3NwZWNpYWwnKTtcbiAgICB9XG4gIH1cbiAgaWYgKGlzVW5kZWZpbmVkKG5hbWUpKSB7XG4gICAgaWYgKGFycmF5ICYmIGtleS5tYXRjaCgvXlxcZCskLykpIHtcbiAgICAgIHJldHVybiBzdHI7XG4gICAgfVxuICAgIG5hbWUgPSBKU09OLnN0cmluZ2lmeSgnJyArIGtleSk7XG4gICAgaWYgKG5hbWUubWF0Y2goL15cIihbYS16QS1aX11bYS16QS1aXzAtOV0qKVwiJC8pKSB7XG4gICAgICBuYW1lID0gbmFtZS5zdWJzdHIoMSwgbmFtZS5sZW5ndGggLSAyKTtcbiAgICAgIG5hbWUgPSBjdHguc3R5bGl6ZShuYW1lLCAnbmFtZScpO1xuICAgIH0gZWxzZSB7XG4gICAgICBuYW1lID0gbmFtZS5yZXBsYWNlKC8nL2csIFwiXFxcXCdcIilcbiAgICAgICAgICAgICAgICAgLnJlcGxhY2UoL1xcXFxcIi9nLCAnXCInKVxuICAgICAgICAgICAgICAgICAucmVwbGFjZSgvKF5cInxcIiQpL2csIFwiJ1wiKTtcbiAgICAgIG5hbWUgPSBjdHguc3R5bGl6ZShuYW1lLCAnc3RyaW5nJyk7XG4gICAgfVxuICB9XG5cbiAgcmV0dXJuIG5hbWUgKyAnOiAnICsgc3RyO1xufVxuXG5cbmZ1bmN0aW9uIHJlZHVjZVRvU2luZ2xlU3RyaW5nKG91dHB1dCwgYmFzZSwgYnJhY2VzKSB7XG4gIHZhciBudW1MaW5lc0VzdCA9IDA7XG4gIHZhciBsZW5ndGggPSBvdXRwdXQucmVkdWNlKGZ1bmN0aW9uKHByZXYsIGN1cikge1xuICAgIG51bUxpbmVzRXN0Kys7XG4gICAgaWYgKGN1ci5pbmRleE9mKCdcXG4nKSA+PSAwKSBudW1MaW5lc0VzdCsrO1xuICAgIHJldHVybiBwcmV2ICsgY3VyLnJlcGxhY2UoL1xcdTAwMWJcXFtcXGRcXGQ/bS9nLCAnJykubGVuZ3RoICsgMTtcbiAgfSwgMCk7XG5cbiAgaWYgKGxlbmd0aCA+IDYwKSB7XG4gICAgcmV0dXJuIGJyYWNlc1swXSArXG4gICAgICAgICAgIChiYXNlID09PSAnJyA/ICcnIDogYmFzZSArICdcXG4gJykgK1xuICAgICAgICAgICAnICcgK1xuICAgICAgICAgICBvdXRwdXQuam9pbignLFxcbiAgJykgK1xuICAgICAgICAgICAnICcgK1xuICAgICAgICAgICBicmFjZXNbMV07XG4gIH1cblxuICByZXR1cm4gYnJhY2VzWzBdICsgYmFzZSArICcgJyArIG91dHB1dC5qb2luKCcsICcpICsgJyAnICsgYnJhY2VzWzFdO1xufVxuXG5cbi8vIE5PVEU6IFRoZXNlIHR5cGUgY2hlY2tpbmcgZnVuY3Rpb25zIGludGVudGlvbmFsbHkgZG9uJ3QgdXNlIGBpbnN0YW5jZW9mYFxuLy8gYmVjYXVzZSBpdCBpcyBmcmFnaWxlIGFuZCBjYW4gYmUgZWFzaWx5IGZha2VkIHdpdGggYE9iamVjdC5jcmVhdGUoKWAuXG5mdW5jdGlvbiBpc0FycmF5KGFyKSB7XG4gIHJldHVybiBBcnJheS5pc0FycmF5KGFyKTtcbn1cbmV4cG9ydHMuaXNBcnJheSA9IGlzQXJyYXk7XG5cbmZ1bmN0aW9uIGlzQm9vbGVhbihhcmcpIHtcbiAgcmV0dXJuIHR5cGVvZiBhcmcgPT09ICdib29sZWFuJztcbn1cbmV4cG9ydHMuaXNCb29sZWFuID0gaXNCb29sZWFuO1xuXG5mdW5jdGlvbiBpc051bGwoYXJnKSB7XG4gIHJldHVybiBhcmcgPT09IG51bGw7XG59XG5leHBvcnRzLmlzTnVsbCA9IGlzTnVsbDtcblxuZnVuY3Rpb24gaXNOdWxsT3JVbmRlZmluZWQoYXJnKSB7XG4gIHJldHVybiBhcmcgPT0gbnVsbDtcbn1cbmV4cG9ydHMuaXNOdWxsT3JVbmRlZmluZWQgPSBpc051bGxPclVuZGVmaW5lZDtcblxuZnVuY3Rpb24gaXNOdW1iZXIoYXJnKSB7XG4gIHJldHVybiB0eXBlb2YgYXJnID09PSAnbnVtYmVyJztcbn1cbmV4cG9ydHMuaXNOdW1iZXIgPSBpc051bWJlcjtcblxuZnVuY3Rpb24gaXNTdHJpbmcoYXJnKSB7XG4gIHJldHVybiB0eXBlb2YgYXJnID09PSAnc3RyaW5nJztcbn1cbmV4cG9ydHMuaXNTdHJpbmcgPSBpc1N0cmluZztcblxuZnVuY3Rpb24gaXNTeW1ib2woYXJnKSB7XG4gIHJldHVybiB0eXBlb2YgYXJnID09PSAnc3ltYm9sJztcbn1cbmV4cG9ydHMuaXNTeW1ib2wgPSBpc1N5bWJvbDtcblxuZnVuY3Rpb24gaXNVbmRlZmluZWQoYXJnKSB7XG4gIHJldHVybiBhcmcgPT09IHZvaWQgMDtcbn1cbmV4cG9ydHMuaXNVbmRlZmluZWQgPSBpc1VuZGVmaW5lZDtcblxuZnVuY3Rpb24gaXNSZWdFeHAocmUpIHtcbiAgcmV0dXJuIGlzT2JqZWN0KHJlKSAmJiBvYmplY3RUb1N0cmluZyhyZSkgPT09ICdbb2JqZWN0IFJlZ0V4cF0nO1xufVxuZXhwb3J0cy5pc1JlZ0V4cCA9IGlzUmVnRXhwO1xuXG5mdW5jdGlvbiBpc09iamVjdChhcmcpIHtcbiAgcmV0dXJuIHR5cGVvZiBhcmcgPT09ICdvYmplY3QnICYmIGFyZyAhPT0gbnVsbDtcbn1cbmV4cG9ydHMuaXNPYmplY3QgPSBpc09iamVjdDtcblxuZnVuY3Rpb24gaXNEYXRlKGQpIHtcbiAgcmV0dXJuIGlzT2JqZWN0KGQpICYmIG9iamVjdFRvU3RyaW5nKGQpID09PSAnW29iamVjdCBEYXRlXSc7XG59XG5leHBvcnRzLmlzRGF0ZSA9IGlzRGF0ZTtcblxuZnVuY3Rpb24gaXNFcnJvcihlKSB7XG4gIHJldHVybiBpc09iamVjdChlKSAmJlxuICAgICAgKG9iamVjdFRvU3RyaW5nKGUpID09PSAnW29iamVjdCBFcnJvcl0nIHx8IGUgaW5zdGFuY2VvZiBFcnJvcik7XG59XG5leHBvcnRzLmlzRXJyb3IgPSBpc0Vycm9yO1xuXG5mdW5jdGlvbiBpc0Z1bmN0aW9uKGFyZykge1xuICByZXR1cm4gdHlwZW9mIGFyZyA9PT0gJ2Z1bmN0aW9uJztcbn1cbmV4cG9ydHMuaXNGdW5jdGlvbiA9IGlzRnVuY3Rpb247XG5cbmZ1bmN0aW9uIGlzUHJpbWl0aXZlKGFyZykge1xuICByZXR1cm4gYXJnID09PSBudWxsIHx8XG4gICAgICAgICB0eXBlb2YgYXJnID09PSAnYm9vbGVhbicgfHxcbiAgICAgICAgIHR5cGVvZiBhcmcgPT09ICdudW1iZXInIHx8XG4gICAgICAgICB0eXBlb2YgYXJnID09PSAnc3RyaW5nJyB8fFxuICAgICAgICAgdHlwZW9mIGFyZyA9PT0gJ3N5bWJvbCcgfHwgIC8vIEVTNiBzeW1ib2xcbiAgICAgICAgIHR5cGVvZiBhcmcgPT09ICd1bmRlZmluZWQnO1xufVxuZXhwb3J0cy5pc1ByaW1pdGl2ZSA9IGlzUHJpbWl0aXZlO1xuXG5leHBvcnRzLmlzQnVmZmVyID0gcmVxdWlyZSgnLi9zdXBwb3J0L2lzQnVmZmVyJyk7XG5cbmZ1bmN0aW9uIG9iamVjdFRvU3RyaW5nKG8pIHtcbiAgcmV0dXJuIE9iamVjdC5wcm90b3R5cGUudG9TdHJpbmcuY2FsbChvKTtcbn1cblxuXG5mdW5jdGlvbiBwYWQobikge1xuICByZXR1cm4gbiA8IDEwID8gJzAnICsgbi50b1N0cmluZygxMCkgOiBuLnRvU3RyaW5nKDEwKTtcbn1cblxuXG52YXIgbW9udGhzID0gWydKYW4nLCAnRmViJywgJ01hcicsICdBcHInLCAnTWF5JywgJ0p1bicsICdKdWwnLCAnQXVnJywgJ1NlcCcsXG4gICAgICAgICAgICAgICdPY3QnLCAnTm92JywgJ0RlYyddO1xuXG4vLyAyNiBGZWIgMTY6MTk6MzRcbmZ1bmN0aW9uIHRpbWVzdGFtcCgpIHtcbiAgdmFyIGQgPSBuZXcgRGF0ZSgpO1xuICB2YXIgdGltZSA9IFtwYWQoZC5nZXRIb3VycygpKSxcbiAgICAgICAgICAgICAgcGFkKGQuZ2V0TWludXRlcygpKSxcbiAgICAgICAgICAgICAgcGFkKGQuZ2V0U2Vjb25kcygpKV0uam9pbignOicpO1xuICByZXR1cm4gW2QuZ2V0RGF0ZSgpLCBtb250aHNbZC5nZXRNb250aCgpXSwgdGltZV0uam9pbignICcpO1xufVxuXG5cbi8vIGxvZyBpcyBqdXN0IGEgdGhpbiB3cmFwcGVyIHRvIGNvbnNvbGUubG9nIHRoYXQgcHJlcGVuZHMgYSB0aW1lc3RhbXBcbmV4cG9ydHMubG9nID0gZnVuY3Rpb24oKSB7XG4gIGNvbnNvbGUubG9nKCclcyAtICVzJywgdGltZXN0YW1wKCksIGV4cG9ydHMuZm9ybWF0LmFwcGx5KGV4cG9ydHMsIGFyZ3VtZW50cykpO1xufTtcblxuXG4vKipcbiAqIEluaGVyaXQgdGhlIHByb3RvdHlwZSBtZXRob2RzIGZyb20gb25lIGNvbnN0cnVjdG9yIGludG8gYW5vdGhlci5cbiAqXG4gKiBUaGUgRnVuY3Rpb24ucHJvdG90eXBlLmluaGVyaXRzIGZyb20gbGFuZy5qcyByZXdyaXR0ZW4gYXMgYSBzdGFuZGFsb25lXG4gKiBmdW5jdGlvbiAobm90IG9uIEZ1bmN0aW9uLnByb3RvdHlwZSkuIE5PVEU6IElmIHRoaXMgZmlsZSBpcyB0byBiZSBsb2FkZWRcbiAqIGR1cmluZyBib290c3RyYXBwaW5nIHRoaXMgZnVuY3Rpb24gbmVlZHMgdG8gYmUgcmV3cml0dGVuIHVzaW5nIHNvbWUgbmF0aXZlXG4gKiBmdW5jdGlvbnMgYXMgcHJvdG90eXBlIHNldHVwIHVzaW5nIG5vcm1hbCBKYXZhU2NyaXB0IGRvZXMgbm90IHdvcmsgYXNcbiAqIGV4cGVjdGVkIGR1cmluZyBib290c3RyYXBwaW5nIChzZWUgbWlycm9yLmpzIGluIHIxMTQ5MDMpLlxuICpcbiAqIEBwYXJhbSB7ZnVuY3Rpb259IGN0b3IgQ29uc3RydWN0b3IgZnVuY3Rpb24gd2hpY2ggbmVlZHMgdG8gaW5oZXJpdCB0aGVcbiAqICAgICBwcm90b3R5cGUuXG4gKiBAcGFyYW0ge2Z1bmN0aW9ufSBzdXBlckN0b3IgQ29uc3RydWN0b3IgZnVuY3Rpb24gdG8gaW5oZXJpdCBwcm90b3R5cGUgZnJvbS5cbiAqL1xuZXhwb3J0cy5pbmhlcml0cyA9IHJlcXVpcmUoJ2luaGVyaXRzJyk7XG5cbmV4cG9ydHMuX2V4dGVuZCA9IGZ1bmN0aW9uKG9yaWdpbiwgYWRkKSB7XG4gIC8vIERvbid0IGRvIGFueXRoaW5nIGlmIGFkZCBpc24ndCBhbiBvYmplY3RcbiAgaWYgKCFhZGQgfHwgIWlzT2JqZWN0KGFkZCkpIHJldHVybiBvcmlnaW47XG5cbiAgdmFyIGtleXMgPSBPYmplY3Qua2V5cyhhZGQpO1xuICB2YXIgaSA9IGtleXMubGVuZ3RoO1xuICB3aGlsZSAoaS0tKSB7XG4gICAgb3JpZ2luW2tleXNbaV1dID0gYWRkW2tleXNbaV1dO1xuICB9XG4gIHJldHVybiBvcmlnaW47XG59O1xuXG5mdW5jdGlvbiBoYXNPd25Qcm9wZXJ0eShvYmosIHByb3ApIHtcbiAgcmV0dXJuIE9iamVjdC5wcm90b3R5cGUuaGFzT3duUHJvcGVydHkuY2FsbChvYmosIHByb3ApO1xufVxuXG59KS5jYWxsKHRoaXMscmVxdWlyZSgnX3Byb2Nlc3MnKSx0eXBlb2YgZ2xvYmFsICE9PSBcInVuZGVmaW5lZFwiID8gZ2xvYmFsIDogdHlwZW9mIHNlbGYgIT09IFwidW5kZWZpbmVkXCIgPyBzZWxmIDogdHlwZW9mIHdpbmRvdyAhPT0gXCJ1bmRlZmluZWRcIiA/IHdpbmRvdyA6IHt9KSJdfQ==
