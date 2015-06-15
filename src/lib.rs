extern crate poisson;
use poisson::{Sample, PoissonDisk};

extern crate nalgebra as na;
pub type Vec2 = na::Vec2<f64>;

extern crate rand;
use rand::Rng;

use std::collections::HashMap;
use std::ops::{Add, Rem};
use std::cmp::{PartialEq, Eq, PartialOrd};
use std::default::Default;
use std::hash::{Hasher, Hash};

#[cfg(test)]
mod test;

pub struct PoissonTileSetGen<R: Rng> {
    poisson: PoissonDisk<R, Vec2>,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct SideCombo(u64, u64, u64, u64);

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct CornerCombo(u64, u64, u64, u64, u64, u64, u64, u64, u64, u64, u64, u64);

pub struct PoissonTileSet {
    tiles: HashMap<CornerCombo, Vec<Sample<Vec2>>>,
    ver_colors: u64,
    hor_colors: u64,
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

impl PoissonTileSet {
    pub fn get_tile<H: StatelessHasher>(&self, hasher: &H, grid_x: i64, grid_y: i64) -> &Vec<Sample<Vec2>> {
        let tl = hasher.hash((grid_x, grid_y + 1));
        let tl_hor = modulo(tl, self.hor_colors);
        let tl_ver = modulo(tl, self.ver_colors);
        let ttlc = modulo(hasher.hash((grid_x, grid_y + 2)), self.hor_colors);
        let tllc = modulo(hasher.hash((grid_x - 1, grid_y + 1)), self.ver_colors);
        let tr = hasher.hash((grid_x + 1, grid_y + 1));
        let tr_hor = modulo(tr, self.hor_colors);
        let tr_ver = modulo(tr, self.ver_colors);
        let ttrc = modulo(hasher.hash((grid_x + 1, grid_y + 2)), self.hor_colors);
        let trrc = modulo(hasher.hash((grid_x + 2, grid_y + 1)), self.ver_colors);
        let bl = hasher.hash((grid_x, grid_y));
        let bl_hor = modulo(bl, self.hor_colors);
        let bl_ver = modulo(bl, self.ver_colors);
        let bblc = modulo(hasher.hash((grid_x, grid_y - 1)), self.hor_colors);
        let bllc = modulo(hasher.hash((grid_x - 1, grid_y)), self.ver_colors);
        let br = hasher.hash((grid_x + 1, grid_y));
        let br_hor = modulo(br, self.hor_colors);
        let br_ver = modulo(br, self.ver_colors);
        let bbrc = modulo(hasher.hash((grid_x + 1, grid_y - 1)), self.hor_colors);
        let brrc = modulo(hasher.hash((grid_x + 2, grid_y)), self.ver_colors);

        let t = modulo(tl_ver + tr_ver, self.ver_colors);
        let r = modulo(tr_hor + br_hor, self.hor_colors);
        let b = modulo(br_ver + bl_ver, self.ver_colors);
        let l = modulo(bl_hor + tl_hor, self.hor_colors);

        let tlt = modulo(tl_hor + ttlc, self.hor_colors);
        let tll = modulo(tl_ver + tllc, self.ver_colors);
        let trt = modulo(tr_hor + ttrc, self.hor_colors);
        let trr = modulo(tr_ver + trrc, self.ver_colors);
        let blb = modulo(bl_hor + bblc, self.hor_colors);
        let bll = modulo(bl_ver + bllc, self.ver_colors);
        let brb = modulo(br_hor + bbrc, self.hor_colors);
        let brr = modulo(br_ver + brrc, self.ver_colors);
        self.tiles.get(&CornerCombo(tll, tlt, trt, trr, brr, brb, blb, bll, t, r, b, l)).unwrap()
    }
}

#[inline]
fn modulo<A, B, C>(n: A, m: B) -> C where A: Rem<B, Output = C>, B: Clone, C: Add<B, Output = C> + Default + PartialOrd {
    let result = n % m.clone();
    if result < C::default() {
        result + m
    } else {
        result
    }
}

impl<R: Rng> PoissonTileSetGen<R> {
    pub fn new(poisson: PoissonDisk<R, Vec2>) -> PoissonTileSetGen<R> {
        PoissonTileSetGen{poisson: poisson}
    }

    pub fn create(&mut self, ver_colors: u64, hor_colors: u64) -> PoissonTileSet {
        println!("Started to create PoissonTileSet...");
        let radius = self.poisson.radius();
        let min = 0.0;//2.0 * radius;
        let scale = 1.0;// + 2.0 * min;
        let max = 1.0 - min;
        let mid = (min + max) / 2.0;
        let side = radius / scale;
        let corner = 2.0 * (side * side / 2.0).sqrt() + radius;
        let hor_edges = self.create_edges(min, mid, max, side, corner, hor_colors, |v| v.pos.y, |v| v.pos.x);
        println!("{} horizontal edges created.", hor_edges.len());
        let ver_edges = self.create_edges(min, mid, max, side, corner, hor_colors, |v| v.pos.x, |v| v.pos.y);
        println!("{} vertical edges created.", ver_edges.len());
        let corners = self.create_corners(min, mid, max, corner, side, hor_colors, &hor_edges, ver_colors, &ver_edges);
        println!("{} corners created.", corners.len());
        let mut tiles = self.create_tiles(min, mid, max, corner, side, hor_colors, &hor_edges, ver_colors, &ver_edges, &corners);
        println!("{} tiles created.", tiles.len());
        for (_, t) in &mut tiles {
            for s in t {
                s.pos = s.pos * scale;
            }
        }
        PoissonTileSet{tiles: tiles, ver_colors: ver_colors, hor_colors: hor_colors}
    }

    fn create_edges<F1, F2>(&mut self, min: f64, mid: f64, max: f64, side: f64, corner: f64, colors: u64, f1: F1, f2: F2) -> HashMap<u64, Vec<Sample<Vec2>>> where F1: Fn(&Sample<Vec2>) -> f64, F2:  Fn(&Sample<Vec2>) -> f64 {
        let mut edges = HashMap::with_capacity(colors as usize);
        for c in 0..colors {
            let mut vectors = vec![];
            self.poisson.create(&mut vectors);
            let v = vectors.clone();
            let edge = vectors.into_iter()
                    .filter(|v| f1(v) > min + corner)
                    .filter(|v| f1(v) < max - corner)
                    .filter(|v| f2(v) < mid + side)
                    .filter(|v| f2(v) > mid - side)
                    .collect();
            debug_edge(c, min, mid, max, side, corner, &v, &edge, &f1);
            edges.insert(c, edge);
        }
        edges
    }

    fn create_corners(&mut self, min: f64, mid: f64, max: f64, corner_size: f64, side_width: f64, hor_colors: u64, hor_edges: &HashMap<u64, Vec<Sample<Vec2>>>, ver_colors: u64, ver_edges: &HashMap<u64, Vec<Sample<Vec2>>>) -> HashMap<SideCombo, Vec<Sample<Vec2>>> {
    	let mut corners = HashMap::with_capacity((hor_colors * hor_colors * ver_colors * ver_colors) as usize);
        for (t, top_set) in hor_edges {
    		let top = new_moved(top_set, Vec2::new(0.0, mid - min));
            for (r, right_set) in ver_edges {
    			let right = new_moved(right_set, Vec2::new(mid - min, 0.0));
                for (b, bot_set) in hor_edges {
    				let bot = new_moved(bot_set, Vec2::new(0.0, -(mid - min)));
    				for (l, left_set) in ver_edges {
    					let left = new_moved(left_set, Vec2::new(-(mid - min), 0.0));
    					self.create_corner(min, mid, max, corner_size, side_width, &mut corners, SideCombo(*t, *r, *b, *l), &top, &bot, &left, &right);
    				}
    			}
    		}
    	}
        corners
    }

    fn create_corner(&mut self, min: f64, mid: f64, max: f64, corner_size: f64, side_width: f64, corners: &mut HashMap<SideCombo, Vec<Sample<Vec2>>>, color: SideCombo, top: &Vec<Sample<Vec2>>, bot: &Vec<Sample<Vec2>>, left: &Vec<Sample<Vec2>>, right: &Vec<Sample<Vec2>>) {
    	let mut combined = Vec::with_capacity(top.len() + bot.len() + left.len() + right.len());
    	combined.extend(top);
    	combined.extend(bot);
    	combined.extend(left);
    	combined.extend(right);
        let mut combined = combined.into_iter().map(|s| *s).collect();
    	self.poisson.create(&mut combined);
        let c = combined.clone();
    	let corner = combined.into_iter()
    			.filter(|v| v.pos.x < mid + corner_size)
    			.filter(|v| v.pos.x > mid - corner_size)
    			.filter(|v| v.pos.y < mid + corner_size)
    			.filter(|v| v.pos.y > mid - corner_size)
    			.filter(|v| v.pos.y < 1.0 * (v.pos.x - (mid - side_width)) + mid + corner_size)
    			.filter(|v| v.pos.y < -1.0 * (v.pos.x - (mid + side_width)) + mid + corner_size)
    			.filter(|v| v.pos.y > 1.0 * (v.pos.x - (mid + side_width)) + mid - corner_size)
    			.filter(|v| v.pos.y > -1.0 * (v.pos.x - (mid - side_width)) + mid - corner_size)
                .collect();
        debug_corner(corners.len(), min, mid, max, corner_size, side_width, color, &top, &bot, &left, &right, &c, &corner);
    	corners.insert(color, corner);
    }


    fn create_tiles(&mut self, min: f64, mid: f64, max: f64, corner_size: f64, side_width: f64, hor_colors: u64, hor_edges: &HashMap<u64, Vec<Sample<Vec2>>>, ver_colors: u64, ver_edges: &HashMap<u64, Vec<Sample<Vec2>>>, corners: &HashMap<SideCombo, Vec<Sample<Vec2>>>) -> HashMap<CornerCombo, Vec<Sample<Vec2>>> {
    	let mut tiles = HashMap::with_capacity((hor_colors.pow(6) + ver_colors.pow(6)) as usize);
        for (t, top_set) in ver_edges {
    		let top = new_moved(top_set, Vec2::new(0.0, mid - min));
            for (r, right_set) in hor_edges {
                let right = new_moved(right_set, Vec2::new(mid - min, 0.0));
                for (b, bot_set) in ver_edges{
                    let bot = new_moved(bot_set, Vec2::new(0.0, -(mid - min)));
    				for (l, left_set) in hor_edges {
                        let left = new_moved(left_set, Vec2::new(-(mid - min), 0.0));
    					self.create_tile_per_corner(min, mid, max, corner_size, side_width, hor_colors, ver_colors, corners, *t, &top, *b, &bot, *l, &left, *r, &right, &mut tiles);
    				}
    			}
    		}
    	}
        tiles
    }

    fn create_tile_per_corner(&mut self, min: f64, mid: f64, max: f64, corner_size: f64, side_width: f64, hor_colors: u64, ver_colors: u64, corners: &HashMap<SideCombo, Vec<Sample<Vec2>>>, t: u64, top: &Vec<Sample<Vec2>>, b: u64, bot: &Vec<Sample<Vec2>>, l: u64, left: &Vec<Sample<Vec2>>, r: u64, right: &Vec<Sample<Vec2>>, tiles: &mut HashMap<CornerCombo, Vec<Sample<Vec2>>>) {
    	for tlt in 0..hor_colors {
            for tll in 0..ver_colors {
                let tlc = SideCombo(tlt, t, l, tll);
    			let tl = corners.get(&tlc).unwrap();
    			let top_left = new_moved(tl, Vec2::new(-(mid - min), mid - min));
                for trt in 0..hor_colors {
                    for trr in 0..ver_colors {
    					let trc = SideCombo(trt, trr, r, t);
    					let tr = corners.get(&trc).unwrap();
    					let top_right = new_moved(tr, Vec2::new(mid - min, mid - min));
    					self.create_tile_per_bottom_corners(min, mid, max, corner_size, side_width, hor_colors, ver_colors, corners, t, top, b, bot, l, left, r, right, tlt, tll, &top_left, trt, trr, &top_right, tiles);
    				}
    			}
    		}
    	}
    }

    fn create_tile_per_bottom_corners(&mut self, min: f64, mid: f64, max: f64, corner_size: f64, side_width: f64, hor_colors: u64, ver_colors: u64, corners: &HashMap<SideCombo, Vec<Sample<Vec2>>>, t: u64, top: &Vec<Sample<Vec2>>, b: u64, bot: &Vec<Sample<Vec2>>, l: u64, left: &Vec<Sample<Vec2>>, r: u64, right: &Vec<Sample<Vec2>>, tlt: u64, tll: u64, top_left: &Vec<Sample<Vec2>>, trt: u64, trr: u64, top_right: &Vec<Sample<Vec2>>, tiles: &mut HashMap<CornerCombo, Vec<Sample<Vec2>>>) {
    	for brb in 0..hor_colors {
            for brr in 0..ver_colors {
    		    let brc = SideCombo(r, brr, brb, b);
    		    let br = corners.get(&brc).unwrap();
    			let bot_right = new_moved(br, Vec2::new(mid - min, -(mid - min)));
                for blb in 0..hor_colors {
                    for bll in 0..ver_colors {
    					let blc = SideCombo(l, b, blb, bll);
    					let bl = corners.get(&blc).unwrap();
    					let bot_left = new_moved(bl, Vec2::new(-(mid - min), -(mid - min)));
    					let cc = CornerCombo(tll, tlt, trt, trr, brr, brb, blb, bll, t, r, b, l);
    					self.create_tile(min, mid, max, corner_size, side_width, cc, top, bot, left, right, top_left, top_right, &bot_right, &bot_left, tiles);
    				}
    			}
    		}
    	}
    }

    fn create_tile(&mut self, min: f64, mid: f64, max: f64, corner_size: f64, side_width: f64, color: CornerCombo, top: &Vec<Sample<Vec2>>, bot: &Vec<Sample<Vec2>>, left: &Vec<Sample<Vec2>>, right: &Vec<Sample<Vec2>>, top_left: &Vec<Sample<Vec2>>, top_right: &Vec<Sample<Vec2>>, bot_right: &Vec<Sample<Vec2>>, bot_left: &Vec<Sample<Vec2>>, tiles: &mut HashMap<CornerCombo, Vec<Sample<Vec2>>>) {
    	let mut combined = Vec::with_capacity(right.len() + left.len() + bot.len() + top.len() + top_left.len() + top_right.len() + bot_right.len() + bot_left.len());
    	combined.extend(right);
    	combined.extend(left);
    	combined.extend(bot);
    	combined.extend(top);
    	combined.extend(top_left);
    	combined.extend(top_right);
    	combined.extend(bot_right);
    	combined.extend(bot_left);
        let mut combined = combined.into_iter().map(|s| *s).collect();
        self.poisson.create(&mut combined);
        let c = combined.clone();
    	let tile = combined.into_iter()
                .filter(|v| v.pos.x >= min)
    			.filter(|v| v.pos.x <= max)
                .filter(|v| v.pos.y >= min)
    			.filter(|v| v.pos.y <= max)
                .collect();
        // debugTile(tiles.len(), min, mid, max, corner_size, side_width, color, &c, &tile, &top, &bot, left, &right, &top_left, &top_right, &bot_right, &bot_left);
    	tiles.insert(color, tile);
    }
}

fn new_moved(vecs: &Vec<Sample<Vec2>>, vec: Vec2) -> Vec<Sample<Vec2>> {
    vecs.iter()
        .map(|s| Sample::new(s.pos + vec, s.radius()))
        .collect()
}

#[cfg(test)]
fn debug_edge<F>(index: u64, min: f64, mid: f64, max: f64, side: f64, corner: f64, vectors: &Vec<Sample<Vec2>>, edge: &Vec<Sample<Vec2>>, f: F) where F: Fn(&Sample<Vec2>) -> f64 {
    if f(&Sample::new(Vec2::new(1.0, 0.0), 0.0)) == 1.0 {
        debug::debug_ver_edge(index, min, mid, max, side, corner, &vectors, &edge);
    } else {
        debug::debug_hor_edge(index, min, mid, max, side, corner, &vectors, &edge);
    }
}
#[cfg(not(test))]
fn debug_edge<F>(index: u64, min: f64, mid: f64, max: f64, side: f64, corner: f64, vectors: &Vec<Sample<Vec2>>, edge: &Vec<Sample<Vec2>>, f: F) where F: Fn(&Sample<Vec2>) -> f64 {}

#[cfg(test)]
fn debug_corner(index: usize, min: f64, mid: f64, max: f64, cornerSize: f64, sideWidth: f64, color: SideCombo, top: &Vec<Sample<Vec2>>, bot: &Vec<Sample<Vec2>>, left: &Vec<Sample<Vec2>>, right: &Vec<Sample<Vec2>>, combined: &Vec<Sample<Vec2>>, corner: &Vec<Sample<Vec2>>){
    debug::debug_corner(index, min, mid, max, cornerSize, sideWidth, color, top, bot, left, right, combined, corner);
}
#[cfg(not(test))]
fn debug_corner(index: usize, min: f64, mid: f64, max: f64, cornerSize: f64, sideWidth: f64, color: SideCombo, top: &Vec<Sample<Vec2>>, bot: &Vec<Sample<Vec2>>, left: &Vec<Sample<Vec2>>, right: &Vec<Sample<Vec2>>, combined: &Vec<Sample<Vec2>>, corner: &Vec<Sample<Vec2>>){}

#[cfg(test)]
fn debugTile(index: usize, min: f64, mid: f64, max: f64, cornerSize: f64, sideWidth: f64, color: CornerCombo, combined: &Vec<Sample<Vec2>>, tile: &Vec<Sample<Vec2>>, top: &Vec<Sample<Vec2>>, bot: &Vec<Sample<Vec2>>, left: &Vec<Sample<Vec2>>, right: &Vec<Sample<Vec2>>, topLeft: &Vec<Sample<Vec2>>, topRight: &Vec<Sample<Vec2>>, botRight: &Vec<Sample<Vec2>>, botLeft: &Vec<Sample<Vec2>>) {
    debug::debugTile(index, min, mid, max, cornerSize, sideWidth, color, combined, tile, top, bot, left, right, topLeft, topRight, botRight, botLeft);
}
#[cfg(not(test))]
fn debugTile(index: usize, min: f64, mid: f64, max: f64, cornerSize: f64, sideWidth: f64, color: CornerCombo, combined: &Vec<Sample<Vec2>>, tile: &Vec<Sample<Vec2>>, top: &Vec<Sample<Vec2>>, bot: &Vec<Sample<Vec2>>, left: &Vec<Sample<Vec2>>, right: &Vec<Sample<Vec2>>, topLeft: &Vec<Sample<Vec2>>, topRight: &Vec<Sample<Vec2>>, botRight: &Vec<Sample<Vec2>>, botLeft: &Vec<Sample<Vec2>>) {}

#[cfg(test)]
mod debug {
    use super::{CornerCombo, SideCombo};

    extern crate image;
    use self::image::{Pixel, ImageBuffer, Rgb};

    use poisson::{Sample, PoissonDisk};

    extern crate nalgebra as na;
    pub type Vec2 = na::Vec2<f64>;
    use na::Norm;

    use std::fs::{self, File};
    use std::path::{Path, PathBuf};
    use std::ops::{DerefMut, Deref};
    use std::collections::HashMap;

    pub fn debug_hor_edge(index: u64, min: f64, mid: f64, max: f64, side: f64, corner: f64, vectors: &Vec<Sample<Vec2>>, edge: &Vec<Sample<Vec2>>) {
        let colors = [
            image::Rgb([0u8, 0u8, 255u8]),
            image::Rgb([255u8, 255u8, 0u8]),
        ];
        let shapes = vec![
        (vec![
            Vec2::new(mid + side, min + corner),
            Vec2::new(mid - side, min + corner),
            Vec2::new(mid - side, max - corner),
            Vec2::new(mid + side, max - corner)], colors[index as usize]),
        // (vec![
        //     Vec2::new(min, min),
        //     Vec2::new(min, max),
        //     Vec2::new(max, max),
        //     Vec2::new(max, min)], image::Rgb([255u8, 255u8, 255u8])),
        ];
        let points = vec![
        (edge.clone(), colors[index as usize]),
        (vectors.iter().filter(|v| !edge.contains(v)).map(|v| v.clone()).collect(), image::Rgb([127u8, 127u8, 127u8]))
        ];
        draw(&*format!("hor{}", index), 256, shapes, points, side);
    }

    pub fn debug_ver_edge(c: u64, min: f64, mid: f64, max: f64, side: f64, corner: f64, vectors: &Vec<Sample<Vec2>>, edge: &Vec<Sample<Vec2>>) {
        let colors = [
            image::Rgb([255u8, 0u8, 0u8]),
            image::Rgb([0u8, 255u8, 0u8]),
        ];
        let shapes = vec![
        (vec![
            Vec2::new(min + corner, mid + side),
            Vec2::new(min + corner, mid - side),
            Vec2::new(max - corner, mid - side),
            Vec2::new(max - corner, mid + side)], colors[c as usize]),
        // (vec![
        //     Vec2::new(min, min),
        //     Vec2::new(min, max),
        //     Vec2::new(max, max),
        //     Vec2::new(max, min)], image::Rgb([255u8, 255u8, 255u8])),
        ];
        let points = vec![
        (edge.clone(), colors[c as usize]),
        (vectors.iter().filter(|v| !edge.contains(v)).map(|v| v.clone()).collect(), image::Rgb([127u8, 127u8, 127u8]))
        ];
        draw(&*format!("ver{}", c), 256, shapes, points, side);
    }

    pub fn debug_corner(index: usize, min: f64, mid: f64, max: f64, cornerSize: f64, sideWidth: f64, color: SideCombo, top: &Vec<Sample<Vec2>>, bot: &Vec<Sample<Vec2>>, left: &Vec<Sample<Vec2>>, right: &Vec<Sample<Vec2>>, combined: &Vec<Sample<Vec2>>, corner: &Vec<Sample<Vec2>>){
        let hor_colors = [
            image::Rgb([0u8, 0u8, 255u8]),
            image::Rgb([255u8, 255u8, 0u8]),
        ];
        let ver_colors = [
            image::Rgb([255u8, 0u8, 0u8]),
            image::Rgb([0u8, 255u8, 0u8]),
        ];
        let lines = vec![
        (vec![
            Vec2::new(mid - cornerSize, mid - sideWidth),
            Vec2::new(mid - cornerSize, mid + sideWidth),
            Vec2::new(mid - sideWidth, mid + cornerSize),
            Vec2::new(mid + sideWidth, mid + cornerSize),
            Vec2::new(mid + cornerSize, mid + sideWidth),
            Vec2::new(mid + cornerSize, mid - sideWidth),
            Vec2::new(mid + sideWidth, mid - cornerSize),
            Vec2::new(mid - sideWidth, mid - cornerSize),
            ], image::Rgb([255u8, 255u8, 255u8])),
        (vec![
            Vec2::new(mid + sideWidth, min + cornerSize + (mid - min)),
            Vec2::new(mid - sideWidth, min + cornerSize + (mid - min)),
            Vec2::new(mid - sideWidth, max - cornerSize + (mid - min)),
            Vec2::new(mid + sideWidth, max - cornerSize + (mid - min))
            ], ver_colors[color.0 as usize]),
        (vec![
            Vec2::new(min + cornerSize + (mid - min), mid + sideWidth),
            Vec2::new(min + cornerSize + (mid - min), mid - sideWidth),
            Vec2::new(max - cornerSize + (mid - min), mid - sideWidth),
            Vec2::new(max - cornerSize + (mid - min), mid + sideWidth)
            ], hor_colors[color.1 as usize]),
        (vec![
            Vec2::new(mid + sideWidth, min + cornerSize - (mid - min)),
            Vec2::new(mid - sideWidth, min + cornerSize - (mid - min)),
            Vec2::new(mid - sideWidth, max - cornerSize - (mid - min)),
            Vec2::new(mid + sideWidth, max - cornerSize - (mid - min))
            ], ver_colors[color.2 as usize]),
        (vec![
            Vec2::new(min + cornerSize - (mid - min), mid + sideWidth),
            Vec2::new(min + cornerSize - (mid - min), mid - sideWidth),
            Vec2::new(max - cornerSize - (mid - min), mid - sideWidth),
            Vec2::new(max - cornerSize - (mid - min), mid + sideWidth)
            ], hor_colors[color.3 as usize]),
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
                    .map(|v| v.clone())
                    .collect(), image::Rgb([127u8, 127u8, 127u8])),
        (top.clone(), ver_colors[color.0 as usize]),
        (bot.clone(), ver_colors[color.2 as usize]),
        (left.clone(), hor_colors[color.3 as usize]),
        (right.clone(), hor_colors[color.1 as usize]),
        ];
		draw(&*format!("corner{}", index), 256, lines, points, sideWidth);
	}

    pub fn debugTile(index: usize, min: f64, mid: f64, max: f64, cornerSize: f64, sideWidth: f64, color: CornerCombo, combined: &Vec<Sample<Vec2>>, tile: &Vec<Sample<Vec2>>, top: &Vec<Sample<Vec2>>, bot: &Vec<Sample<Vec2>>, left: &Vec<Sample<Vec2>>, right: &Vec<Sample<Vec2>>, topLeft: &Vec<Sample<Vec2>>, topRight: &Vec<Sample<Vec2>>, botRight: &Vec<Sample<Vec2>>, botLeft: &Vec<Sample<Vec2>>) {
        let hor_colors = [
            image::Rgb([0u8, 0u8, 255u8]),
            image::Rgb([255u8, 255u8, 0u8]),
        ];
        let ver_colors = [
            image::Rgb([255u8, 0u8, 0u8]),
            image::Rgb([0u8, 255u8, 0u8]),
        ];
        let points = vec![
        (tile.iter()
                .filter(|v| !top.contains(v))
                .filter(|v| !right.contains(v))
                .filter(|v| !bot.contains(v))
                .filter(|v| !left.contains(v))
                .filter(|v| !topLeft.contains(v))
                .filter(|v| !topRight.contains(v))
                .filter(|v| !botRight.contains(v))
                .filter(|v| !botLeft.contains(v))
                .map(|v| v.clone())
                .collect(), image::Rgb([255u8, 255u8, 255u8])),
        (top.iter()
                .filter(|v| tile.contains(v))
                .map(|v| v.clone())
                .collect(), ver_colors[color.8 as usize]),
        (right.iter()
                .filter(|v| tile.contains(v))
                .map(|v| v.clone())
                .collect(), hor_colors[color.9 as usize]),
        (bot.iter()
                .filter(|v| tile.contains(v))
                .map(|v| v.clone())
                .collect(), ver_colors[color.10 as usize]),
        (left.iter()
                .filter(|v| tile.contains(v))
                .map(|v| v.clone())
                .collect(), hor_colors[color.11 as usize]),
        (top.iter()
                .filter(|v| !tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[ver_colors[color.8 as usize], image::Rgb([127u8, 127u8, 127u8])])),
        (right.iter()
                .filter(|v| !tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[hor_colors[color.9 as usize], image::Rgb([127u8, 127u8, 127u8])])),
        (bot.iter()
                .filter(|v| !tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[ver_colors[color.10 as usize], image::Rgb([127u8, 127u8, 127u8])])),
        (left.iter()
                .filter(|v| !tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[hor_colors[color.11 as usize], image::Rgb([127u8, 127u8, 127u8])])),
        (topLeft.iter()
                .filter(|v| tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[hor_colors[color.0 as usize], ver_colors[color.1 as usize]])),
        (topRight.iter()
                .filter(|v| tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[hor_colors[color.3 as usize], ver_colors[color.2 as usize]])),
        (botRight.iter()
                .filter(|v| tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[hor_colors[color.4 as usize], ver_colors[color.5 as usize]])),
        (botLeft.iter()
                .filter(|v| tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[hor_colors[color.7 as usize], ver_colors[color.6 as usize]])),
        (topLeft.iter()
                .filter(|v| !tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[combine(&[hor_colors[color.0 as usize], ver_colors[color.1 as usize]]), image::Rgb([127u8, 127u8, 127u8])])),
        (topRight.iter()
                .filter(|v| !tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[combine(&[hor_colors[color.3 as usize], ver_colors[color.2 as usize]]), image::Rgb([127u8, 127u8, 127u8])])),
        (botRight.iter()
                .filter(|v| !tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[combine(&[hor_colors[color.4 as usize], ver_colors[color.5 as usize]]), image::Rgb([127u8, 127u8, 127u8])])),
        (botLeft.iter()
                .filter(|v| !tile.contains(v))
                .map(|v| v.clone())
                .collect(), combine(&[combine(&[hor_colors[color.7 as usize], ver_colors[color.6 as usize]]), image::Rgb([127u8, 127u8, 127u8])])),
        (combined.iter()
                .filter(|v| !tile.contains(v))
                .filter(|v| !top.contains(v))
                .filter(|v| !right.contains(v))
                .filter(|v| !bot.contains(v))
                .filter(|v| !left.contains(v))
                .filter(|v| !topLeft.contains(v))
                .filter(|v| !topRight.contains(v))
                .filter(|v| !botRight.contains(v))
                .filter(|v| !botLeft.contains(v))
                .map(|v| v.clone())
                .collect(), image::Rgb([127u8, 127u8, 127u8]))
        ];
        // (tll, tlt, trt, trr, brr, brb, blb, bll, t, r, b, l)
        let lines = vec![
        (vec![
            Vec2::new(min + cornerSize, max),
            Vec2::new(min + cornerSize, max - sideWidth),
            Vec2::new(max - cornerSize, max - sideWidth),
            Vec2::new(max - cornerSize, max)
            ], ver_colors[color.0 as usize]),
        (vec![
            Vec2::new(max, min + cornerSize),
            Vec2::new(max - sideWidth, min + cornerSize),
            Vec2::new(max - sideWidth, max - cornerSize),
            Vec2::new(max, max - cornerSize)
            ], hor_colors[color.1 as usize]),
        (vec![
            Vec2::new(min + cornerSize, min + sideWidth),
            Vec2::new(min + cornerSize, min),
            Vec2::new(max - cornerSize, min),
            Vec2::new(max - cornerSize, min + sideWidth)
            ], ver_colors[color.2 as usize]),
        (vec![
            Vec2::new(min + sideWidth, min + cornerSize),
            Vec2::new(min, min + cornerSize),
            Vec2::new(min, max - cornerSize),
            Vec2::new(min + sideWidth, max - cornerSize)
            ], hor_colors[color.3 as usize]),
        (vec![
            Vec2::new(min, min),
            Vec2::new(min, min + cornerSize),
            Vec2::new(min + sideWidth, min + cornerSize),
            Vec2::new(min + cornerSize, min + sideWidth),
            Vec2::new(min + cornerSize, min)
            ], combine(&[hor_colors[color.7 as usize], ver_colors[color.6 as usize]])),
        (vec![
            Vec2::new(min, max),
            Vec2::new(min + cornerSize, max),
            Vec2::new(min + cornerSize, max - sideWidth),
            Vec2::new(min + sideWidth, max - cornerSize),
            Vec2::new(min, max - cornerSize)
            ], combine(&[hor_colors[color.0 as usize], ver_colors[color.1 as usize]])),
        (vec![
            Vec2::new(max - cornerSize, max - sideWidth),
            Vec2::new(max - cornerSize, max),
            Vec2::new(max, max),
            Vec2::new(max, max - cornerSize),
            Vec2::new(max - sideWidth, max - cornerSize)
            ], combine(&[hor_colors[color.3 as usize], ver_colors[color.2 as usize]])),
        (vec![
            Vec2::new(max - cornerSize, min),
            Vec2::new(max - cornerSize, min + sideWidth),
            Vec2::new(max - sideWidth, min + cornerSize),
            Vec2::new(max, min + cornerSize),
            Vec2::new(max, min)
            ], combine(&[hor_colors[color.4 as usize], ver_colors[color.5 as usize]])),
        ];
        draw(&*format!("tile{}", index), 256, lines, points, sideWidth);
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

    fn is_legal_poisson(vecs: & Vec<(Vec<Sample<Vec2>>, Rgb<u8>)>, vec: &Sample<Vec2>) -> bool {
        for &(ref vv, _) in vecs {
            for v in vv {
                if vec.pos == v.pos {
                    continue;
                }
                let diff = vec.pos - v.pos;
                let dist = diff.norm();
                let allowed_dist = v.radius() + vec.radius();
                if dist < allowed_dist {
                    return false;
                }
            }
        }
        true
    }

    fn draw(path: &str, side: u32, shapes: Vec<(Vec<Vec2>, Rgb<u8>)>, points: Vec<(Vec<Sample<Vec2>>, Rgb<u8>)>, radius: f64) {
        let mut imgbuf = image::ImageBuffer::new(side, side);
    	let r = (side as f64 * radius) as i32;
    	let r2 = r * r;
        for &(ref s, ref c) in &points {
            for v in s {
                let x = (v.pos.x * side as f64) as i32;
    			let y = (v.pos.y * side as f64) as i32;
                let legal = is_legal_poisson(&points, &v);
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
    			let mut lastX = (vector.x * side as f64) as i32;
    			let mut lastY = (vector.y * side as f64) as i32;
                for v in shape {
                    let x = (v.x * side as f64) as i32;
    				let y = (v.y * side as f64) as i32;
    				line(&mut imgbuf, side, lastX, lastY, x, y, color);
    				lastX = x;
    				lastY = y;
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
