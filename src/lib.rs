#[macro_use]
extern crate lazy_static;

extern crate poisson;
use poisson::{Sample, PoissonGen};

extern crate nalgebra as na;
use na::{Vec2, Pnt2, FloatPnt};

extern crate rand;
use rand::Rng;

use std::collections::HashMap;
use std::ops::{Add, Rem};
use std::cmp::{Eq, PartialOrd};
use std::default::Default;
use std::hash::{Hasher, Hash};
use std::borrow::{Borrow, Cow};

#[cfg(test)]
mod test;

pub type VecSample = Vec<Sample<Vec2<f64>>>;

pub struct PoissonTileSetGen<R: Rng> {
    poisson: PoissonGen<R, Vec2<f64>>,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct SideCombo {
    t: u64,
    r: u64,
    b: u64,
    l: u64,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct TileColor {
    tll: u64,
    tlt: u64,
    trt: u64,
    trr: u64,
    brr: u64,
    brb: u64,
    blb: u64,
    bll: u64,
    sides: SideCombo,
}

pub trait StatelessHasher {
    fn hash<H: Hash>(&self, target: H) -> u64;
}

impl<T: Hasher + Clone> StatelessHasher for T {
    fn hash<H: Hash>(&self, target: H) -> u64{
        let mut hasher = self.clone();
        target.hash(&mut hasher);
        hasher.finish()
    }
}

pub trait Mod<RHS=Self> {
    type Output;
    fn modulo(self, rhs: RHS) -> Self::Output;
}

impl<A, B, C> Mod<B> for A where A: Rem<B, Output=C>, B: Clone, C: Add<B, Output=C> + Default + PartialOrd {
    type Output = C;
    fn modulo(self, rhs: B) -> Self::Output {
        let result = self % rhs.clone();
        if result < Self::Output::default() {
            result + rhs
        } else {
            result
        }
    }
}

pub struct Tile<'a> {
    color: TileColor,
    samples: Cow<'a, VecSample>,
    coords: Pnt2<i64>,
}

impl<'a> Tile<'a> {
    fn find_closest(&self, pnt: Pnt2<f64>) -> (Vec2<f64>, f64, usize) {
        let mut index = 0;
        let mut dist = std::f64::MAX;
        let mut sample = None;
        for (i, s) in self.samples.iter().enumerate() {
            let d = pnt.dist(s.pos.as_pnt());
            if d < dist {
                dist = d;
                index = i;
                sample = Some(s.pos);
            }
        }
        (sample.unwrap(), dist, index)
    }
}

pub struct PoissonTileSet {
    tiles: HashMap<TileColor, VecSample>,
    radius: f64,
    ver_colors: u64,
    hor_colors: u64,
}

impl PoissonTileSet {
    pub fn get_tile<H: StatelessHasher>(&self, hasher: &H, grid_coord: Pnt2<i64>) -> Tile {
        let mut tile = self.get_nonmoved_tile(hasher, grid_coord);
        tile.samples = Cow::Owned(new_moved(tile.samples.borrow(), Vec2::new(grid_coord.x as f64, grid_coord.y as f64)));
        tile
    }

    pub fn get_nonmoved_tile<H: StatelessHasher>(&self, hasher: &H, grid_coord: Pnt2<i64>) -> Tile {
        let hor = self.hor_colors;
        let ver = self.ver_colors;

        let bl = hasher.hash(grid_coord + Vec2::new(0, 0));
        let bl_hor = bl.modulo(hor);
        let bl_ver = bl.modulo(ver);
        let bblc = hasher.hash(grid_coord + Vec2::new(0, -1)).modulo(hor);
        let bllc = hasher.hash(grid_coord + Vec2::new(-1, 0)).modulo(ver);

        let tl = hasher.hash(grid_coord + Vec2::new(0, 1));
        let tl_hor = tl.modulo(hor);
        let tl_ver = tl.modulo(ver);
        let ttlc = hasher.hash(grid_coord + Vec2::new(0, 2)).modulo(hor);
        let tllc = hasher.hash(grid_coord + Vec2::new(-1, 1)).modulo(ver);

        let br = hasher.hash(grid_coord + Vec2::new(1, 0));
        let br_hor = br.modulo(hor);
        let br_ver = br.modulo(ver);
        let bbrc = hasher.hash(grid_coord + Vec2::new(1, -1)).modulo(hor);
        let brrc = hasher.hash(grid_coord + Vec2::new(2, 0)).modulo(ver);

        let tr = hasher.hash(grid_coord + Vec2::new(1, 1));
        let tr_hor = tr.modulo(hor);
        let tr_ver = tr.modulo(ver);
        let ttrc = hasher.hash(grid_coord + na::Vec2::new(1, 2)).modulo(hor);
        let trrc = hasher.hash(grid_coord + na::Vec2::new(2, 1)).modulo(ver);

        let color = TileColor{
            tll: (tl_ver + tllc).modulo(ver),
            tlt: (tl_hor + ttlc).modulo(hor),
            trr: (tr_ver + trrc).modulo(ver),
            trt: (tr_hor + ttrc).modulo(hor),
            brr: (br_ver + brrc).modulo(ver),
            brb: (br_hor + bbrc).modulo(hor),
            bll: (bl_ver + bllc).modulo(ver),
            blb: (bl_hor + bblc).modulo(hor),
            sides: SideCombo {
                t: (tl_ver + tr_ver).modulo(ver),
                r: (tr_hor + br_hor).modulo(hor),
                b: (br_ver + bl_ver).modulo(ver),
                l: (bl_hor + tl_hor).modulo(hor),
            },
        };
        Tile{
            samples: Cow::Borrowed(self.tiles.get(&color).unwrap()),
            color: color,
            coords: grid_coord,
        }
    }

    //TODO: This breaks on tile corners...
    pub fn find_closest<H: StatelessHasher>(&self, hasher: &H, pnt: Pnt2<f64>) -> (Vec2<f64>, u64, f64, bool) {
        let grid_coord = Pnt2::new(pnt.x.floor(), pnt.y.floor());
        let grid = Pnt2::new(grid_coord.x as i64, grid_coord.y as i64);
        let mut tiles = vec![self.get_tile(hasher, grid)];
        let y_diff = pnt.y - grid_coord.y;
        let x_diff = pnt.x - grid_coord.x;
        // radius * 2 smallest number that is larger than largest empty circle possible
        let largest_empty_circle = self.radius * 2.;
        let largest_empty_circle_sqr = largest_empty_circle.powi(2);

        // tiles.push(self.get_tile(hasher, grid + Vec2::new(-1, -1)));
        // tiles.push(self.get_tile(hasher, grid + Vec2::new(-1, 0)));
        // tiles.push(self.get_tile(hasher, grid + Vec2::new(-1, 1)));
        // tiles.push(self.get_tile(hasher, grid + Vec2::new(0, -1)));
        //
        // tiles.push(self.get_tile(hasher, grid + Vec2::new(0, 1)));
        // tiles.push(self.get_tile(hasher, grid + Vec2::new(1, -1)));
        // tiles.push(self.get_tile(hasher, grid + Vec2::new(1, 0)));
        // tiles.push(self.get_tile(hasher, grid + Vec2::new(1, 1)));

        if pnt.sqdist(&(grid_coord + Vec2::new(0., 0.))) < largest_empty_circle_sqr {
            tiles.push(self.get_tile(hasher, grid + Vec2::new(-1, 0)));
            tiles.push(self.get_tile(hasher, grid + Vec2::new(0, -1)));
            tiles.push(self.get_tile(hasher, grid + Vec2::new(-1, -1)));
            // left, bot, letf-bot
        } else if pnt.sqdist(&(grid_coord + Vec2::new(0., 1.))) < largest_empty_circle_sqr {
            tiles.push(self.get_tile(hasher, grid + Vec2::new(-1, 0)));
            tiles.push(self.get_tile(hasher, grid + Vec2::new(0, 1)));
            tiles.push(self.get_tile(hasher, grid + Vec2::new(-1, 1)));
            // left, top, left-top
        } else if pnt.sqdist(&(grid_coord + Vec2::new(1., 0.))) < largest_empty_circle_sqr {
            tiles.push(self.get_tile(hasher, grid + Vec2::new(1, 0)));
            tiles.push(self.get_tile(hasher, grid + Vec2::new(0, -1)));
            tiles.push(self.get_tile(hasher, grid + Vec2::new(1, -1)));
            // right, bot, rigth-bot
        } else if pnt.sqdist(&(grid_coord + Vec2::new(1., 1.))) < largest_empty_circle_sqr {
            tiles.push(self.get_tile(hasher, grid + Vec2::new(1, 0)));
            tiles.push(self.get_tile(hasher, grid + Vec2::new(0, 1)));
            tiles.push(self.get_tile(hasher, grid + Vec2::new(1, 1)));
            // right, top, right-top
        } else if y_diff < largest_empty_circle {
            tiles.push(self.get_tile(hasher, grid + Vec2::new(0, -1)));
            // bot, this
        } else if y_diff > 1. - largest_empty_circle {
            tiles.push(self.get_tile(hasher, grid + Vec2::new(0, 1)));
            // top, this
        } else if x_diff < largest_empty_circle {
            tiles.push(self.get_tile(hasher, grid + Vec2::new(-1, 0)));
            // left, this
        } else if x_diff > 1. - largest_empty_circle {
            tiles.push(self.get_tile(hasher, grid + Vec2::new(1, 0)));
            // right, this
        }

        let mut closest_dist = std::f64::MAX;
        let mut closest = None;
        let mut id = 0;
        let mut inside = false;
        for tile in tiles {
            let (cur, dist, index) = tile.find_closest(pnt);
            if dist < closest_dist {
                closest = Some(cur);
                closest_dist = dist;
                id = hasher.hash((tile.coords, index));
                if dist <= self.radius {
                    inside = true;
                    break;
                }
            }
        }
        (closest.unwrap(), id, closest_dist, inside)
    }
}

#[derive(Copy, Clone)]
pub struct Dimensions {
    min: f64,
    mid: f64,
    max: f64,
    side: f64,
    corner: f64,
}

impl Dimensions {
    fn new(min: f64, max: f64, radius: f64) -> Self {
        let max = max - min;
        let mid = (min + max) / 2.;
        let corner = (2. * radius.powi(2)).sqrt() + radius;//2.0 * (side.powi(2) / 2.0).sqrt() + radius;
        Dimensions{min: min, mid: mid, max: max, side: radius, corner: corner}
    }
}

impl<R: Rng> PoissonTileSetGen<R> {
    pub fn new(poisson: PoissonGen<R, Vec2<f64>>) -> PoissonTileSetGen<R> {
        PoissonTileSetGen{poisson: poisson}
    }

    pub fn create(&mut self, ver_colors: u64, hor_colors: u64) -> PoissonTileSet {
        let dims = Dimensions::new(0.0, 1.0, self.poisson.radius());
        println!("Started to create PoissonTileSet...");

        let hor_edges = self.create_edges(dims, hor_colors, |v| v.y, |v| v.x);
        println!("{} horizontal edges created.", hor_edges.len());

        let ver_edges = self.create_edges(dims, hor_colors, |v| v.x, |v| v.y);
        println!("{} vertical edges created.", ver_edges.len());

        let mut corners = HashMap::with_capacity((hor_colors.pow(2) * ver_colors.pow(2)) as usize);
        Self::for_each_side_color(dims, &hor_edges, &ver_edges, |c, t, b, l, r| self.create_corner(dims, &mut corners, c, t, b, l, r));
        println!("{} corners created.", corners.len());

        let mut tiles = HashMap::with_capacity((hor_colors.pow(6) + ver_colors.pow(6)) as usize);
        Self::for_each_side_color(dims, &ver_edges, &hor_edges, |c, t, b, l, r| self.create_tile_per_corner(dims, hor_colors, ver_colors, &corners, &mut tiles, c, t, b, l, r));
        println!("{} tiles created.", tiles.len());

        PoissonTileSet {
            tiles: tiles,
            ver_colors: ver_colors,
            hor_colors: hor_colors,
            radius: self.poisson.radius(),
        }
    }

    fn create_edges<F1, F2>(&mut self, dims: Dimensions, colors: u64, f1: F1, f2: F2) -> HashMap<u64, VecSample> where F1: Fn(&Vec2<f64>) -> f64, F2:  Fn(&Vec2<f64>) -> f64 {
        let mut edges = HashMap::with_capacity(colors as usize);
        for c in 0..colors {
            let mut edge = vec![];
            self.poisson.generate_with_limit(&mut edge, |v| {
                    f1(&v) > dims.min + dims.corner &&
                    f1(&v) < dims.max - dims.corner &&
                    f2(&v) < dims.mid + dims.side &&
                    f2(&v) > dims.mid - dims.side
                });
            debug_edge(c, dims, &edge, &edge, &f1);
            edges.insert(c, edge);
        }
        edges
    }

    fn for_each_side_color<F>(dims: Dimensions, edges1: &HashMap<u64, VecSample>, edges2: &HashMap<u64, VecSample>, mut f: F) where F: FnMut(SideCombo, &VecSample, &VecSample, &VecSample, &VecSample) {
        let half = dims.mid - dims.min;
        for (t, top) in edges1 {
    		let top = new_moved(top, Vec2::new(0., half));
            for (r, right) in edges2 {
    			let right = new_moved(right, Vec2::new(half, 0.));
                for (b, bot) in edges1 {
    				let bot = new_moved(bot, Vec2::new(0., -half));
    				for (l, left) in edges2 {
    					let left = new_moved(left, Vec2::new(-half, 0.));
                        (f)(SideCombo{t: *t, r: *r, b: *b, l: *l}, &top, &bot, &left, &right);
    				}
    			}
    		}
    	}
    }

    fn create_corner(&mut self, dims: Dimensions, corners: &mut HashMap<SideCombo, VecSample>, color: SideCombo, top: &VecSample, bot: &VecSample, left: &VecSample, right: &VecSample) {
    	let mut combined = Vec::with_capacity(top.len() + bot.len() + left.len() + right.len());
    	combined.extend(top);
    	combined.extend(bot);
    	combined.extend(left);
    	combined.extend(right);
        let mut combined = combined.into_iter().cloned().collect();
    	self.poisson.generate_with_limit(&mut combined, |v| {
                v.x < dims.mid + dims.corner &&
                v.x > dims.mid - dims.corner &&
                v.y < dims.mid + dims.corner &&
                v.y > dims.mid - dims.corner &&
                v.y < 1. * (v.x - (dims.mid - dims.side)) + dims.mid + dims.corner &&
                v.y < -1. * (v.x - (dims.mid + dims.side)) + dims.mid + dims.corner &&
                v.y > 1. * (v.x - (dims.mid + dims.side)) + dims.mid - dims.corner &&
                v.y > -1. * (v.x - (dims.mid - dims.side)) + dims.mid - dims.corner
            });
        let c = combined.clone();
    	let corner = combined.into_iter()
                            .filter(|v| !top.contains(v))
                            .filter(|v| !bot.contains(v))
                            .filter(|v| !left.contains(v))
                            .filter(|v| !right.contains(v))
                            .collect();
        debug_corner(corners.len(), dims, color, &top, &bot, &left, &right, &c, &corner);
    	corners.insert(color, corner);
    }

    fn create_tile_per_corner(&mut self, dims: Dimensions, hor_colors: u64, ver_colors: u64, corners: &HashMap<SideCombo, VecSample>, tiles: &mut HashMap<TileColor, VecSample>, col: SideCombo, top: &VecSample, bot: &VecSample, left: &VecSample, right: &VecSample) {
    	let half = dims.mid - dims.min;
        for tlt in 0..hor_colors {
            for tll in 0..ver_colors {
                let tlc = SideCombo{t: tlt, r: col.t, b: col.l, l: tll};
    			let tl = corners.get(&tlc).unwrap();
    			let top_left = new_moved(tl, Vec2::new(-half, half));
                for trt in 0..hor_colors {
                    for trr in 0..ver_colors {
    					let trc = SideCombo{t: trt, r: trr, b: col.r, l: col.t};
    					let tr = corners.get(&trc).unwrap();
    					let top_right = new_moved(tr, Vec2::new(half, half));
    					self.create_tile_per_bottom_corners(dims, hor_colors, ver_colors, corners, col, top, bot, left, right, tlt, tll, &top_left, trt, trr, &top_right, tiles);
    				}
    			}
    		}
    	}
    }

    fn create_tile_per_bottom_corners(&mut self, dims: Dimensions, hor_colors: u64, ver_colors: u64, corners: &HashMap<SideCombo, VecSample>, col: SideCombo, top: &VecSample, bot: &VecSample, left: &VecSample, right: &VecSample, tlt: u64, tll: u64, top_left: &VecSample, trt: u64, trr: u64, top_right: &VecSample, tiles: &mut HashMap<TileColor, VecSample>) {
    	let half = dims.mid - dims.min;
        for brb in 0..hor_colors {
            for brr in 0..ver_colors {
    		    let brc = SideCombo{t: col.r, r: brr, b: brb, l: col.b};
    		    let br = corners.get(&brc).unwrap();
    			let bot_right = new_moved(br, Vec2::new(half, -half));
                for blb in 0..hor_colors {
                    for bll in 0..ver_colors {
    					let blc = SideCombo{t: col.l, r: col.b, b: blb, l: bll};
    					let bl = corners.get(&blc).unwrap();
    					let bot_left = new_moved(bl, Vec2::new(-half, -half));
    					let cc = TileColor{tll: tll, tlt: tlt, trt: trt, trr: trr, brr: brr, brb: brb, blb: blb, bll: bll, sides: col};//t: col.t, r: col.r, b: col.b, l: col.l};
    					self.create_tile(dims, cc, top, bot, left, right, top_left, top_right, &bot_right, &bot_left, tiles);
    				}
    			}
    		}
    	}
    }

    fn create_tile(&mut self, dims: Dimensions, color: TileColor, top: &VecSample, bot: &VecSample, left: &VecSample, right: &VecSample, top_left: &VecSample, top_right: &VecSample, bot_right: &VecSample, bot_left: &VecSample, tiles: &mut HashMap<TileColor, VecSample>) {
        let mut combined = Vec::with_capacity(right.len() + left.len() + bot.len() + top.len() + top_left.len() + top_right.len() + bot_right.len() + bot_left.len());
    	combined.extend(right);
    	combined.extend(left);
    	combined.extend(bot);
    	combined.extend(top);
    	combined.extend(top_left);
    	combined.extend(top_right);
    	combined.extend(bot_right);
    	combined.extend(bot_left);
        let mut combined = combined.into_iter().cloned().collect();
        self.poisson.generate_with_limit(&mut combined, |v| {
                Some(v).iter()
                    // .filter(|v| v.x >= dims.min + dims.side)
                    // .filter(|v| v.x <= dims.max - dims.side)
                    // .filter(|v| v.y >= dims.min + dims.side)
                    // .filter(|v| v.y <= dims.max - dims.side)
                    // .filter(|v| v.y < 1f64 * (v.x - (dims.mid - dims.side)) + dims.mid + dims.corner)
                    // .filter(|v| v.y < -1f64 * (v.x - (dims.mid + dims.side)) + dims.mid + dims.corner)
                    // .filter(|v| v.y > 1f64 * (v.x - (dims.mid + dims.side)) + dims.mid - dims.corner)
                    // .filter(|v| v.y > -1f64 * (v.x - (dims.mid - dims.side)) + dims.mid - dims.corner)
                    .next().is_some()
            });
        let c = combined.clone();
        let mut tile = combined.into_iter()
                                .filter(|v| v.pos.x >= dims.min)
                                .filter(|v| v.pos.x <= dims.max)
                                .filter(|v| v.pos.y >= dims.min)
                                .filter(|v| v.pos.y <= dims.max)
                                .collect::<VecSample>();
        tile.sort_by(|v, o| {
            match v.pos.x.partial_cmp(&o.pos.x).unwrap() {
                std::cmp::Ordering::Equal => {
                    v.pos.y.partial_cmp(&o.pos.y).unwrap()
                }
                r @ _ => {
                    r
                }
            }
        });

        // let scale = (dims.max - dims.min - dims.side) / (dims.max - dims.min);
        // let mut combined: Vec<Sample<Vec2>> = combined.into_iter()
        //     .map(|s: &Sample<Vec2>| Sample::new(s.pos * scale, s.radius() * scale))
        //     .cloned()
        //     .collect();
        // let cc: Vec<Sample<Vec2>> = combined.clone();
        // //TODO: Scaling helps, but there still can be overlaps in corners and it leaves holes into the distribution
        // let radius = self.poisson.radius();
        // self.poisson.set_radius(radius * scale);
        // self.poisson.generate(&mut combined);
        // self.poisson.set_radius(radius);
        // let c = combined.clone();
    	// let mut tile: Vec<Sample<Vec2>> = combined.into_iter()
        //         .map(|s| (cc.contains(&s), Sample::new(s.pos / scale, s.radius() / scale)))
        //         .filter(|v| {
        //             if v.0 {
        //                 v.1.pos.x >= dims.min &&
        //                 v.1.pos.x <= dims.max &&
        //                 v.1.pos.y >= dims.min &&
        //                 v.1.pos.y <= dims.max
        //             } else {
        //                 v.1.pos.x >= dims.min + dims.side &&
        //                 v.1.pos.x <= dims.max - dims.side &&
        //                 v.1.pos.y >= dims.min + dims.side &&
        //                 v.1.pos.y <= dims.max - dims.side
        //             }})
        //         .map(|s| s.1)
        //         // .filter(|v| v.pos.y < 1.0 * (v.pos.x - (mid - side_width)) + mid + corner_size)
        //         // .filter(|v| v.pos.y < -1.0 * (v.pos.x - (mid + side_width)) + mid + corner_size)
        //         // .filter(|v| v.pos.y > 1.0 * (v.pos.x - (mid + side_width)) + mid - corner_size)
        //         // .filter(|v| v.pos.y > -1.0 * (v.pos.x - (mid - side_width)) + mid - corner_size)
        //         .collect();
        // tile.extend(right);
        // tile.extend(left);
        // tile.extend(bot);
        // tile.extend(top);
        // tile.extend(top_left);
        // tile.extend(top_right);
        // tile.extend(bot_right);
        // tile.extend(bot_left);
        // let mut tile = combined.into_iter()
        //         .filter(|v| v.pos.x >= min)
        //         .filter(|v| v.pos.x <= max)
        //         .filter(|v| v.pos.y >= min)
        //         .filter(|v| v.pos.y <= max)
        //         .collect();

        // TODO: debug_tile(tiles.len(), dims, color, &c, &tile, &top, &bot, left, &right, &top_left, &top_right, &bot_right, &bot_left);
    	tiles.insert(color, tile);
    }
}

pub fn new_moved(vecs: &VecSample, vec: Vec2<f64>) -> VecSample {
    vecs.iter()
        .map(|s| Sample::new(s.pos + vec, s.radius()))
        .collect()
}

#[cfg(test)]
fn debug_edge<F>(index: u64, dims: Dimensions, vectors: &VecSample, edge: &VecSample, f: F) where F: Fn(&Vec2<f64>) -> f64 {
    if f(&Vec2::new(1., 0.)) == 1. {
        debug::debug_ver_edge(index, dims, &vectors, &edge);
    } else {
        debug::debug_hor_edge(index, dims, &vectors, &edge);
    }
}
#[cfg(not(test))]
#[allow(unused_variables)]
fn debug_edge<F>(index: u64, dims: Dimensions, vectors: &VecSample, edge: &VecSample, f: F) where F: Fn(&Vec2<f64>) -> f64 {}

#[cfg(test)]
fn debug_corner(index: usize, dims: Dimensions, color: SideCombo, top: &VecSample, bot: &VecSample, left: &VecSample, right: &VecSample, combined: &VecSample, corner: &VecSample){
    debug::debug_corner(index, dims, color, top, bot, left, right, combined, corner);
}
#[cfg(not(test))]
#[allow(unused_variables)]
fn debug_corner(index: usize, dims: Dimensions, color: SideCombo, top: &VecSample, bot: &VecSample, left: &VecSample, right: &VecSample, combined: &VecSample, corner: &VecSample){}

#[cfg(test)]
fn debug_tile(index: usize, dims: Dimensions, color: TileColor, combined: &VecSample, tile: &VecSample, top: &VecSample, bot: &VecSample, left: &VecSample, right: &VecSample, top_left: &VecSample, top_right: &VecSample, bot_right: &VecSample, bot_left: &VecSample) {
    debug::debug_tile(index, dims, color, combined, tile, top, bot, left, right, top_left, top_right, bot_right, bot_left);
}
#[cfg(not(test))]
#[allow(unused_variables)]
fn debug_tile(index: usize, dims: Dimensions, color: TileColor, combined: &VecSample, tile: &VecSample, top: &VecSample, bot: &VecSample, left: &VecSample, right: &VecSample, top_left: &VecSample, top_right: &VecSample, bot_right: &VecSample, bot_left: &VecSample) {}

#[cfg(test)]
mod debug {
    use super::*;

    extern crate image;
    use self::image::{ImageBuffer, Rgb};

    extern crate rand;

    pub use poisson::{Sample, PoissonDisk};

    extern crate nalgebra as na;
    use na::{Vec2, Norm, Pnt2};

    use std::fs::{self, File};
    use std::path::PathBuf;
    use std::ops::{DerefMut, Deref};
    use std::hash::SipHasher;
    use std::borrow::{Borrow, Cow};

    lazy_static! {
        static ref HOR_COLORS: [Rgb<u8>; 2] = [
            image::Rgb([0u8, 0u8, 255u8]),
            image::Rgb([255u8, 255u8, 0u8]),
        ];
        static ref VER_COLORS: [Rgb<u8>; 2] = [
            image::Rgb([255u8, 0u8, 0u8]),
            image::Rgb([0u8, 255u8, 0u8]),
        ];
    }

    pub fn debug_hor_edge(index: u64, d: Dimensions, vectors: &VecSample, edge: &VecSample) {
        let shapes = vec![
        (vec![
            Vec2::new(d.mid + d.side, d.min + d.corner),
            Vec2::new(d.mid - d.side, d.min + d.corner),
            Vec2::new(d.mid - d.side, d.max - d.corner),
            Vec2::new(d.mid + d.side, d.max - d.corner)], HOR_COLORS[index as usize]),
        // (vec![
        //     Vec2::new(min, min),
        //     Vec2::new(min, max),
        //     Vec2::new(max, max),
        //     Vec2::new(max, min)], image::Rgb([255u8, 255u8, 255u8])),
        ];
        let points = vec![
        (edge.clone(), HOR_COLORS[index as usize]),
        (vectors.iter().filter(|v| !edge.contains(v)).cloned().collect(), image::Rgb([127u8, 127u8, 127u8]))
        ];
        draw(&*format!("hor{}", index), |_, _| {}, 256, shapes, points, d.side);
    }

    pub fn debug_ver_edge(c: u64, d: Dimensions, vectors: &VecSample, edge: &VecSample) {
        let shapes = vec![
        (vec![
            Vec2::new(d.min + d.corner, d.mid + d.side),
            Vec2::new(d.min + d.corner, d.mid - d.side),
            Vec2::new(d.max - d.corner, d.mid - d.side),
            Vec2::new(d.max - d.corner, d.mid + d.side)], VER_COLORS[c as usize]),
        // (vec![
        //     Vec2::new(min, min),
        //     Vec2::new(min, max),
        //     Vec2::new(max, max),
        //     Vec2::new(max, min)], image::Rgb([255u8, 255u8, 255u8])),
        ];
        let points = vec![
        (edge.clone(), VER_COLORS[c as usize]),
        (vectors.iter().filter(|v| !edge.contains(v)).cloned().collect(), image::Rgb([127u8, 127u8, 127u8]))
        ];
        draw(&*format!("ver{}", c), |_, _| {}, 256, shapes, points, d.side);
    }

    pub fn debug_corner(index: usize, d: Dimensions, color: SideCombo, top: &VecSample, bot: &VecSample, left: &VecSample, right: &VecSample, combined: &VecSample, corner: &VecSample){
        let lines = vec![
        (vec![
            Vec2::new(d.mid - d.corner, d.mid - d.side),
            Vec2::new(d.mid - d.corner, d.mid + d.side),
            Vec2::new(d.mid - d.side, d.mid + d.corner),
            Vec2::new(d.mid + d.side, d.mid + d.corner),
            Vec2::new(d.mid + d.corner, d.mid + d.side),
            Vec2::new(d.mid + d.corner, d.mid - d.side),
            Vec2::new(d.mid + d.side, d.mid - d.corner),
            Vec2::new(d.mid - d.side, d.mid - d.corner),
            ], image::Rgb([255u8, 255u8, 255u8])),
        (vec![
            Vec2::new(d.mid + d.side, d.min + d.corner + (d.mid - d.min)),
            Vec2::new(d.mid - d.side, d.min + d.corner + (d.mid - d.min)),
            Vec2::new(d.mid - d.side, d.max - d.corner + (d.mid - d.min)),
            Vec2::new(d.mid + d.side, d.max - d.corner + (d.mid - d.min))
            ], HOR_COLORS[color.t as usize]),
        (vec![
            Vec2::new(d.min + d.corner + (d.mid - d.min), d.mid + d.side),
            Vec2::new(d.min + d.corner + (d.mid - d.min), d.mid - d.side),
            Vec2::new(d.max - d.corner + (d.mid - d.min), d.mid - d.side),
            Vec2::new(d.max - d.corner + (d.mid - d.min), d.mid + d.side)
            ], VER_COLORS[color.r as usize]),
        (vec![
            Vec2::new(d.mid + d.side, d.min + d.corner - (d.mid - d.min)),
            Vec2::new(d.mid - d.side, d.min + d.corner - (d.mid - d.min)),
            Vec2::new(d.mid - d.side, d.max - d.corner - (d.mid - d.min)),
            Vec2::new(d.mid + d.side, d.max - d.corner - (d.mid - d.min))
            ], HOR_COLORS[color.b as usize]),
        (vec![
            Vec2::new(d.min + d.corner - (d.mid - d.min), d.mid + d.side),
            Vec2::new(d.min + d.corner - (d.mid - d.min), d.mid - d.side),
            Vec2::new(d.max - d.corner - (d.mid - d.min), d.mid - d.side),
            Vec2::new(d.max - d.corner - (d.mid - d.min), d.mid + d.side)
            ], VER_COLORS[color.l as usize]),
        // (vec![
        //     Vec2::new(min, min),
        //     Vec2::new(min, max),
        //     Vec2::new(max, max),
        //     Vec2::new(max, min)], image::Rgb([255u8, 255u8, 255u8])),
        ];
        let points = vec![
        (corner.clone(), image::Rgb([255u8, 255u8, 255u8])),
        (combined.iter()
                    .filter(|v| !corner.contains(v))
                    .filter(|v| !top.contains(v))
                    .filter(|v| !right.contains(v))
                    .filter(|v| !bot.contains(v))
                    .filter(|v| !left.contains(v))
                    .cloned()
                    .collect(), image::Rgb([127u8, 127u8, 127u8])),
        (top.clone(), HOR_COLORS[color.t as usize]),
        (bot.clone(), HOR_COLORS[color.b as usize]),
        (left.clone(), VER_COLORS[color.l as usize]),
        (right.clone(), VER_COLORS[color.r as usize]),
        ];
		draw(&*format!("corner{}", index), |_, _| {}, 256, lines, points, d.side);
	}

    pub fn debug_tile(index: usize, d: Dimensions, color: TileColor, combined: &VecSample, tile: &VecSample, top: &VecSample, bot: &VecSample, left: &VecSample, right: &VecSample, top_left: &VecSample, top_right: &VecSample, bot_right: &VecSample, bot_left: &VecSample) {
        let points = vec![
        (tile.iter()
                .filter(|v| !top.contains(v))
                .filter(|v| !right.contains(v))
                .filter(|v| !bot.contains(v))
                .filter(|v| !left.contains(v))
                .filter(|v| !top_left.contains(v))
                .filter(|v| !top_right.contains(v))
                .filter(|v| !bot_right.contains(v))
                .filter(|v| !bot_left.contains(v))
                .cloned()
                .collect(), image::Rgb([255u8, 255u8, 255u8])),
        (top.iter()
                .filter(|v| tile.contains(v))
                .cloned()
                .collect(), VER_COLORS[color.sides.t as usize]),
        (right.iter()
                .filter(|v| tile.contains(v))
                .cloned()
                .collect(), HOR_COLORS[color.sides.r as usize]),
        (bot.iter()
                .filter(|v| tile.contains(v))
                .cloned()
                .collect(), VER_COLORS[color.sides.b as usize]),
        (left.iter()
                .filter(|v| tile.contains(v))
                .cloned()
                .collect(), HOR_COLORS[color.sides.l as usize]),
        (top.iter()
                .filter(|v| !tile.contains(v))
                .cloned()
                .collect(), combine(&[VER_COLORS[color.sides.t as usize], image::Rgb([127u8, 127u8, 127u8])])),
        (right.iter()
                .filter(|v| !tile.contains(v))
                .cloned()
                .collect(), combine(&[HOR_COLORS[color.sides.r as usize], image::Rgb([127u8, 127u8, 127u8])])),
        (bot.iter()
                .filter(|v| !tile.contains(v))
                .cloned()
                .collect(), combine(&[VER_COLORS[color.sides.b as usize], image::Rgb([127u8, 127u8, 127u8])])),
        (left.iter()
                .filter(|v| !tile.contains(v))
                .cloned()
                .collect(), combine(&[HOR_COLORS[color.sides.l as usize], image::Rgb([127u8, 127u8, 127u8])])),
        (top_left.iter()
                .filter(|v| tile.contains(v))
                .cloned()
                .collect(), combine(&[HOR_COLORS[color.tll as usize], VER_COLORS[color.tlt as usize]])),
        (top_right.iter()
                .filter(|v| tile.contains(v))
                .cloned()
                .collect(), combine(&[HOR_COLORS[color.trt as usize], VER_COLORS[color.trr as usize]])),
        (bot_right.iter()
                .filter(|v| tile.contains(v))
                .cloned()
                .collect(), combine(&[HOR_COLORS[color.brr as usize], VER_COLORS[color.brb as usize]])),
        (bot_left.iter()
                .filter(|v| tile.contains(v))
                .cloned()
                .collect(), combine(&[HOR_COLORS[color.blb as usize], VER_COLORS[color.bll as usize]])),
        (top_left.iter()
                .filter(|v| !tile.contains(v))
                .cloned()
                .collect(), combine(&[combine(&[HOR_COLORS[color.tll as usize], VER_COLORS[color.tlt as usize]]), image::Rgb([127u8, 127u8, 127u8])])),
        (top_right.iter()
                .filter(|v| !tile.contains(v))
                .cloned()
                .collect(), combine(&[combine(&[HOR_COLORS[color.trt as usize], VER_COLORS[color.trr as usize]]), image::Rgb([127u8, 127u8, 127u8])])),
        (bot_right.iter()
                .filter(|v| !tile.contains(v))
                .cloned()
                .collect(), combine(&[combine(&[HOR_COLORS[color.brr as usize], VER_COLORS[color.brb as usize]]), image::Rgb([127u8, 127u8, 127u8])])),
        (bot_left.iter()
                .filter(|v| !tile.contains(v))
                .cloned()
                .collect(), combine(&[combine(&[HOR_COLORS[color.blb as usize], VER_COLORS[color.bll as usize]]), image::Rgb([127u8, 127u8, 127u8])])),
        (combined.iter()
                .filter(|v| !tile.contains(v))
                .filter(|v| !top.contains(v))
                .filter(|v| !right.contains(v))
                .filter(|v| !bot.contains(v))
                .filter(|v| !left.contains(v))
                .filter(|v| !top_left.contains(v))
                .filter(|v| !top_right.contains(v))
                .filter(|v| !bot_right.contains(v))
                .filter(|v| !bot_left.contains(v))
                .cloned()
                .collect(), image::Rgb([127u8, 127u8, 127u8]))
        ];
        // (tll, tlt, trt, trr, brr, brb, blb, bll, t, r, b, l)
        let lines = vec![
        (vec![
            Vec2::new(d.min + d.corner, d.max),
            Vec2::new(d.min + d.corner, d.max - d.side),
            Vec2::new(d.max - d.corner, d.max - d.side),
            Vec2::new(d.max - d.corner, d.max)
            ], VER_COLORS[color.tll as usize]),
        (vec![
            Vec2::new(d.max, d.min + d.corner),
            Vec2::new(d.max - d.side, d.min + d.corner),
            Vec2::new(d.max - d.side, d.max - d.corner),
            Vec2::new(d.max, d.max - d.corner)
            ], HOR_COLORS[color.tlt as usize]),
        (vec![
            Vec2::new(d.min + d.corner, d.min + d.side),
            Vec2::new(d.min + d.corner, d.min),
            Vec2::new(d.max - d.corner, d.min),
            Vec2::new(d.max - d.corner, d.min + d.side)
            ], VER_COLORS[color.trt as usize]),
        (vec![
            Vec2::new(d.min + d.side, d.min + d.corner),
            Vec2::new(d.min, d.min + d.corner),
            Vec2::new(d.min, d.max - d.corner),
            Vec2::new(d.min + d.side, d.max - d.corner)
            ], HOR_COLORS[color.trr as usize]),
        (vec![
            Vec2::new(d.min, d.min),
            Vec2::new(d.min, d.min + d.corner),
            Vec2::new(d.min + d.side, d.min + d.corner),
            Vec2::new(d.min + d.corner, d.min + d.side),
            Vec2::new(d.min + d.corner, d.min)
            ], combine(&[HOR_COLORS[color.blb as usize], VER_COLORS[color.brb as usize]])),
        (vec![
            Vec2::new(d.min, d.max),
            Vec2::new(d.min + d.corner, d.max),
            Vec2::new(d.min + d.corner, d.max - d.side),
            Vec2::new(d.min + d.side, d.max - d.corner),
            Vec2::new(d.min, d.max - d.corner)
            ], combine(&[HOR_COLORS[color.tll as usize], VER_COLORS[color.tlt as usize]])),
        (vec![
            Vec2::new(d.max - d.corner, d.max - d.side),
            Vec2::new(d.max - d.corner, d.max),
            Vec2::new(d.max, d.max),
            Vec2::new(d.max, d.max - d.corner),
            Vec2::new(d.max - d.side, d.max - d.corner)
            ], combine(&[HOR_COLORS[color.trr as usize], VER_COLORS[color.trt as usize]])),
        (vec![
            Vec2::new(d.max - d.corner, d.min),
            Vec2::new(d.max - d.corner, d.min + d.side),
            Vec2::new(d.max - d.side, d.min + d.corner),
            Vec2::new(d.max, d.min + d.corner),
            Vec2::new(d.max, d.min)
            ], combine(&[HOR_COLORS[color.brr as usize], VER_COLORS[color.brb as usize]])),
        ];
        draw(&*format!("tile{}", index),  |_, _| {}, 256, lines, points, d.side);
    }

    #[test]
    fn test_if_it_works() {
        let poisson = PoissonDisk::new(rand::weak_rng()).build_samples::<Vec2<f64>>(25, 0.8);
        let side_width = poisson.radius();
        let corner_size = 1.5 * side_width;
        let mut gen = PoissonTileSetGen::new(poisson);
        let tile_set = gen.create(2, 2);
        let hasher = SipHasher::new();
        let mut samples = vec![];
        let size = 4;
        let mut lines = vec![];
        for x in -size..size {
            for y in -size..size {
                let m = Vec2::new((size + x) as f64, (size + y) as f64);
                let Tile{color: c, samples: tile, ..} = tile_set.get_tile(&hasher, Pnt2::new(x, y));
                let mut ss = new_moved(tile.borrow(), Vec2::new(size as f64, size as f64));
                for s in &mut ss {
                    s.pos = s.pos / (size * 2) as f64;
                }
                let col = combine(&[VER_COLORS[c.tll as usize], image::Rgb([255u8, 255u8, 255u8])]);
                let col = combine(&[HOR_COLORS[c.tlt as usize], col, col]);
                let col = combine(&[HOR_COLORS[c.trt as usize], col, col]);
                let col = combine(&[HOR_COLORS[c.trr as usize], col, col]);
                let col = combine(&[HOR_COLORS[c.brr as usize], col, col]);
                let col = combine(&[HOR_COLORS[c.brb as usize], col, col]);
                let col = combine(&[HOR_COLORS[c.blb as usize], col, col]);
                let col = combine(&[HOR_COLORS[c.bll as usize], col, col]);
                let col = combine(&[HOR_COLORS[c.sides.t as usize], col, col]);
                let col = combine(&[HOR_COLORS[c.sides.r as usize], col, col]);
                let col = combine(&[HOR_COLORS[c.sides.b as usize], col, col]);
                let col = combine(&[HOR_COLORS[c.sides.l as usize], col, col]);
                samples.push((ss, col));
                lines.push((vec![
                    Vec2::new(-corner_size, -side_width) + m,
                    Vec2::new(-corner_size, side_width) + m,
                    Vec2::new(-side_width, corner_size) + m,
                    Vec2::new(side_width, corner_size) + m,
                    Vec2::new(corner_size, side_width) + m,
                    Vec2::new(corner_size, -side_width) + m,
                    Vec2::new(side_width, -corner_size) + m,
                    Vec2::new(-side_width, -corner_size) + m
                    ], image::Rgb([255u8, 255u8, 255u8])));
                lines.push((vec![
                    Vec2::new(corner_size, 1.0) + m,
                    Vec2::new(corner_size, 1.0 - side_width) + m,
                    Vec2::new(1.0 - corner_size, 1.0 - side_width) + m,
                    Vec2::new(1.0 - corner_size, 1.0) + m
                    ], VER_COLORS[c.sides.t as usize]));
                lines.push((vec![
                    Vec2::new(1.0, corner_size) + m,
                    Vec2::new(1.0 - side_width, corner_size) + m,
                    Vec2::new(1.0 - side_width, 1.0 - corner_size) + m,
                    Vec2::new(1.0, 1.0 - corner_size) + m
                    ], HOR_COLORS[c.sides.r as usize]));
                lines.push((vec![
                    Vec2::new(corner_size, 0.0 + side_width) + m,
                    Vec2::new(corner_size, 0.0) + m,
                    Vec2::new(1.0 - corner_size, 0.0) + m,
                    Vec2::new(1.0 - corner_size, 0.0 + side_width) + m
                    ], VER_COLORS[c.sides.b as usize]));
                lines.push((vec![
                    Vec2::new(0.0 + side_width, corner_size) + m,
                    Vec2::new(0.0, corner_size) + m,
                    Vec2::new(0.0, 1.0 - corner_size) + m,
                    Vec2::new(0.0 + side_width, 1.0 - corner_size) + m
                    ], HOR_COLORS[c.sides.l as usize]));
            }
        }
        let lines = lines.into_iter()
             .map(|(l, c)| (l.into_iter().map(|v| v / (size * 2) as f64).collect(), c))
             .collect();
        // for &mut (ref mut l, _) in &mut lines {
        //     for ref mut v in l {
        //         v = **v / (size * 2) as f64;
        //     }
        // }
        let side = size as u32 * 2 * 256;
        let mut imgbuf = image::ImageBuffer::new(side, side);
        for x in 0..side {
            for y in 0..side {
                //TODO: Coordinates are wrong??? Or is the find_closest broken?
                let xx = x as f64 / (side as f64 / (size * 2) as f64) - size as f64;
                let yy = y as f64 / (side as f64 / (size * 2) as f64) - size as f64;
                let (sample, id, dist, inside) = tile_set.find_closest(&hasher, Pnt2::new(xx, yy));
                let mut col = image::Rgb([(id >> 16 & 0xFF) as u8, (id >> 8  & 0xFF)as u8, (id & 0xFF) as u8]);
                if inside {
                    col = image::Rgb([0, 0, 0]);
                }
                pixel(&mut imgbuf, side, x as i32, y as i32, col);
            }
        }
        let mut p = PathBuf::new();
        p.push("visualise");
        let _ = fs::create_dir(p.clone());
        p.push("tiling_func");
        let ref mut fout = File::create(p.with_extension("png")).unwrap();
        let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);

        draw(&*format!("tiling"), |imgbuf, side| {
            for x in 0..side {
                for y in 0..side {
                    //TODO: Coordinates are wrong??? Or is the find_closest broken?
                    let xx = x as f64 / (side as f64 / (size * 2) as f64) - size as f64;
                    let yy = y as f64 / (side as f64 / (size * 2) as f64) - size as f64;
                    let (sample, id, dist, inside) = tile_set.find_closest(&hasher, Pnt2::new(xx, yy));
                    let mut col = image::Rgb([(id >> 16 & 0xFF) as u8, (id >> 8  & 0xFF)as u8, (id & 0xFF) as u8]);
                    if inside {
                        col = image::Rgb([0, 0, 0]);
                    }
                    pixel(imgbuf, side, x as i32, y as i32, col);
                }
            }
        }, size as u32 * 2 * 256, lines, samples, side_width / (size * 2) as f64);
        // visualise(&samples, size as u32 * 2, 64, "tiling");
        // println!("{:?}", tile_set.get_tile(&hasher, 0, 0));
    }

    fn combine(colors: &[Rgb<u8>]) -> Rgb<u8> {
        let mut red = 0;
        let mut green = 0;
        let mut blue = 0;
        for color in colors {
            red += color.data[0] as u64;
            green += color.data[1] as u64;
            blue += color.data[2] as u64;
        }
        image::Rgb([(red / colors.len() as u64) as u8, (green / colors.len() as u64) as u8, (blue / colors.len() as u64) as u8])
    }

    fn is_legal_poisson(vecs: &Vec<(VecSample, Rgb<u8>)>, vec: &Sample<Vec2<f64>>, radius: f64) -> bool {
        for &(ref vv, _) in vecs {
            for v in vv {
                if vec.pos == v.pos {
                    continue;
                }
                let diff = vec.pos - v.pos;
                let dist = diff.norm();
                let allowed_dist = radius * 2.;//v.radius() + vec.radius();
                if dist < allowed_dist {
                    return false;
                }
            }
        }
        true
    }

    fn draw<F>(path: &str, mut f: F, side: u32, shapes: Vec<(Vec<Vec2<f64>>, Rgb<u8>)>, points: Vec<(VecSample, Rgb<u8>)>, radius: f64) where F: FnMut(&mut ImageBuffer<Rgb<u8>, Vec<u8>>, u32) {
        let mut imgbuf = image::ImageBuffer::new(side, side);
        (f)(&mut imgbuf, side);
    	let r = (side as f64 * radius) as i32;
    	let r2 = r * r;
        for &(ref s, ref c) in &points {
            for v in s {
                let x = (v.pos.x * side as f64) as i32;
    			let y = (v.pos.y * side as f64) as i32;
                let legal = is_legal_poisson(&points, &v, radius);
                for xx in -r..(r + 1) {
                    for yy in -r..(r + 1) {
                        let squared = xx * xx + yy * yy;
                        if squared < r2 && (legal || squared > r2 / 2) {
                            pixel(&mut imgbuf, side, (x + xx) as i32, (y + yy) as i32, *c);
                        }
                    }
                }
            }
        }

        for (shape, color) in shapes {
            if !shape.is_empty() {
    			let vector = shape[shape.len() - 1];
    			let mut last_x = (vector.x * side as f64) as i32;
    			let mut last_y = (vector.y * side as f64) as i32;
                for v in shape {
                    let x = (v.x * side as f64) as i32;
    				let y = (v.y * side as f64) as i32;
    				line(&mut imgbuf, side, last_x, last_y, x, y, color);
    				last_x = x;
    				last_y = y;
                }
    		}
        }
        let mut p = PathBuf::new();
        p.push("visualise");
        let _ = fs::create_dir(p.clone());
        p.push(path);
        let ref mut fout = File::create(p.with_extension("png")).unwrap();
        let _ = image::ImageRgb8(imgbuf).save(fout, image::PNG);
    }

    fn line<C>(imgbuf: &mut ImageBuffer<Rgb<u8>, C>, side: u32, x0: i32, y0: i32, x1: i32, y1: i32, color: Rgb<u8>) where C: Deref<Target=[u8]> + DerefMut {
        if x1 == x0 {
            if y0 < y1 {
                for y in y0..(y1 + 1) {
                    pixel(imgbuf, side, x1, y, color);
                }
            } else {
                for y in y1..(y0 + 1) {
                    pixel(imgbuf, side, x1, y, color);
                }
            }
            return;
        }
        let deltax = (x1 - x0) as f64;
        let deltay  = (y1 - y0) as f64;
        let mut error = 0.0;
        let deltaerr = (deltay / deltax).abs();
        let mut y = y0;
        let sign = if y0 < y1 { 1 } else { -1 };
        if x0 < x1 {
            for x in x0..(x1 + 1) {
                pixel(imgbuf, side, x, y, color);
                error += deltaerr;
                while error >= 0.5 {
                    pixel(imgbuf, side, x, y, color);
                    y += sign;
                    error -= 1.0;
                }
            }
        } else {
            for x in (x1..(x0 + 1)).rev() {
                pixel(imgbuf, side, x, y, color);
                error += deltaerr;
                while error >= 0.5 {
                    pixel(imgbuf, side, x, y, color);
                    y += sign;
                    error -= 1.0;
                }
            }
        }

    }

    fn pixel<C>(imgbuf: &mut ImageBuffer<Rgb<u8>, C>, side: u32, x: i32, y: i32, color: Rgb<u8>) where C: Deref<Target=[u8]> + DerefMut {
        if x >= side as i32 || x < 0 {
            return;
        }
        if y >= side as i32 || y < 0 {
            return;
        }
        imgbuf.put_pixel(x as u32, y as u32, color);
    }
}
