extern crate poisson;
pub use poisson::{Sample, PoissonDisk};

extern crate nalgebra as na;
pub type Vec2 = na::Vec2<f64>;

extern crate rand;

extern crate image;

use ::PoissonTileSetGen;

use std::hash::SipHasher;
use std::fs::{self, File};
use std::path::{Path, PathBuf};

#[test]
fn test_if_it_works() {
    let poisson = PoissonDisk::<_, Vec2>::with_samples(rand::weak_rng(), 25, 0.8, false);
    let mut gen = PoissonTileSetGen::new(poisson);
    let tile_set = gen.create(2, 2);
    let hasher = SipHasher::new();
    let mut samples = vec![];
    let size = 8;
    for x in -size..size {
        for y in -size..size {
            let tile = tile_set.get_tile(&hasher, x, y);
            samples.extend(new_moved(tile, Vec2::new((size + x) as f64, (size + y) as f64)));
        }
    }
    for s in &mut samples {
        s.pos = s.pos / (size * 2) as f64;
    }
    visualise(&samples, size as u32 * 2, 64, "tiling");
    // println!("{:?}", tile_set.get_tile(&hasher, 0, 0));
}

fn new_moved(vecs: &Vec<Sample<Vec2>>, vec: Vec2) -> Vec<Sample<Vec2>> {
    vecs.iter()
        .map(|s| Sample::new(s.pos + vec, s.radius()))
        .collect()
}

fn visualise(samples: &Vec<Sample<Vec2>>, tiles: u32, size_per_tile: u32, path: &str) {
    let size = tiles * size_per_tile;
    let mut imgbuf = image::ImageBuffer::new(size, size);
    let color = image::Rgb([255 as u8, 255 as u8, 255 as u8]);
    let middle = image::Rgb([255 as u8, 0 as u8, 0 as u8]);
    for sample in samples {
        let xx = (sample.pos.x * size as f64) as i32;
        let yy = (sample.pos.y * size as f64) as i32;
        let radius = (sample.radius() * size_per_tile as f64) as i32;
        for x in -radius..(radius + 1) {
            for y in -radius..(radius + 1) {
                if x * x + y * y < radius * radius {
                    let color = if x == 0 || y == 0 {
                        middle
                    } else {
                        color
                    };
                    // let cur_x = if periodicity {
                    //     modulo(xx + x, size as i32)
                    // } else {
                        let cur_x = xx + x;
                        if cur_x < 0 || cur_x >= size as i32 {
                            continue;
                        }
                        // cur_x
                    // };
                    // let cur_y = if periodicity {
                    //     modulo(yy + y, size as i32)
                    // } else {
                        let cur_y = yy + y;
                        if cur_y < 0 || cur_y >= size as i32 {
                            continue;
                        }
                    //     cur_y
                    // };
                    imgbuf.put_pixel(cur_x as u32, cur_y as u32, color);
                }
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
