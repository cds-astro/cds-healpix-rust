#![cfg_attr(test, feature(test))]

#[cfg(test)]
extern crate test;

#[macro_use]
extern crate serde_derive;

extern crate cdshealpix;

use std::io::prelude::*;
use std::fs::File;
use std::time::{Duration, Instant};

use serde_derive::{Deserialize, Serialize};
use serde::{Deserialize, Serialize};
use serde_json::Result;

use cdshealpix::nested::bmoc::*;

#[derive(Serialize, Deserialize)]
pub struct Sdss {
  #[serde(rename = "1")]
  d1: Vec<u64>,
  #[serde(rename = "2")]
  d2: Vec<u64>,
  #[serde(rename = "3")]
  d3: Vec<u64>,
  #[serde(rename = "4")]
  d4: Vec<u64>,
  #[serde(rename = "5")]
  d5: Vec<u64>,
  #[serde(rename = "6")]
  d6: Vec<u64>,
  #[serde(rename = "7")]
  d7: Vec<u64>,
  #[serde(rename = "8")]
  d8: Vec<u64>,
  #[serde(rename = "9")]
  d9: Vec<u64>,
  #[serde(rename = "10")]
  d10: Vec<u64>,
  #[serde(rename = "11")]
  d11: Vec<u64>,
  #[serde(rename = "12")]
  d12: Vec<u64>,
}

impl Sdss {
  
  pub fn to_bmoc(&self) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(12, self.d1.len() + self.d2.len() + self.d3.len() + self.d4.len() 
      + self.d5.len() + self.d6.len() + self.d7.len() + self.d8.len() + self.d9.len() 
      + self.d10.len() + self.d11.len() + self.d12.len());
    for h in &self.d1 {
      builder.push(1, *h, true);
    }
    for h in &self.d2 {
      builder.push(2, *h, true);
    }
    for h in &self.d3 {
      builder.push(3, *h, true);
    }
    for h in &self.d4 {
      builder.push(4, *h, true);
    }
    for h in &self.d5 {
      builder.push(5, *h, true);
    }
    for h in &self.d6 {
      builder.push(6, *h, true);
    }
    for h in &self.d7 {
      builder.push(7, *h, true);
    }
    for h in &self.d8 {
      builder.push(8, *h, true);
    }
    for h in &self.d9 {
      builder.push(9, *h, true);
    }
    for h in &self.d10 {
      builder.push(10, *h, true);
    }
    for h in &self.d11 {
      builder.push(11, *h, true);
    }
    for h in &self.d12 {
      builder.push(12, *h, true);
    }
    builder.to_bmoc_from_unordered()
  }
  
}

pub fn load_sdss() -> Result<Sdss> {
  let mut file = File::open("resources/sdss.moc.json").expect("Unable to open the SDSS file");
  let mut contents = String::new();
  file.read_to_string(&mut contents).expect("Unable to read the SDSS file");
  let sdss: Sdss = serde_json::from_str(&contents)?;
  Ok(sdss)
}

pub fn load_sdss_not() -> Result<Sdss> {
  let mut file = File::open("resources/sdss_not.moc.json").expect("Unable to open the SDSS not file");
  let mut contents = String::new();
  file.read_to_string(&mut contents).expect("Unable to read the SDSS not file");
  let sdss: Sdss = serde_json::from_str(&contents)?;
  Ok(sdss)
}

pub fn load_s_or_g() -> Result<Sdss> {
  let mut file = File::open("resources/sdss_or_glimpse.moc.json").expect("Unable to open the SDSS OR GLIMPSE not file");
  let mut contents = String::new();
  file.read_to_string(&mut contents).expect("Unable to read the SDSS OR GLIMPSE not file");
  let sdss: Sdss = serde_json::from_str(&contents)?;
  Ok(sdss)
}

pub fn load_s_xor_g() -> Result<Sdss> {
  let mut file = File::open("resources/sdss_xor_glimpse.moc.json").expect("Unable to open the SDSS XOR GLIMPSE not file");
  let mut contents = String::new();
  file.read_to_string(&mut contents).expect("Unable to read the SDSS XOR GLIMPSE not file");
  let sdss: Sdss = serde_json::from_str(&contents)?;
  Ok(sdss)
}


#[derive(Serialize, Deserialize)]
pub struct Glimpse {
  #[serde(rename = "4")]
  d4: Vec<u64>,
  #[serde(rename = "5")]
  d5: Vec<u64>,
  #[serde(rename = "6")]
  d6: Vec<u64>,
  #[serde(rename = "7")]
  d7: Vec<u64>,
  #[serde(rename = "8")]
  d8: Vec<u64>,
  #[serde(rename = "9")]
  d9: Vec<u64>,
}

impl Glimpse {

  pub fn to_bmoc(&self) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(9, self.d4.len() + self.d5.len() + self.d6.len() 
      + self.d7.len() + self.d8.len() + self.d9.len());
    for h in &self.d4 {
      builder.push(4, *h, true);
    }
    for h in &self.d5 {
      builder.push(5, *h, true);
    }
    for h in &self.d6 {
      builder.push(6, *h, true);
    }
    for h in &self.d7 {
      builder.push(7, *h, true);
    }
    for h in &self.d8 {
      builder.push(8, *h, true);
    }
    for h in &self.d9 {
      builder.push(9, *h, true);
    }
    builder.to_bmoc_from_unordered()
  }

}

pub fn load_glimpse() -> Result<Glimpse> {
  let mut file = File::open("resources/glimpse.moc.json").expect("Unable to open the Glimpse file");
  let mut contents = String::new();
  file.read_to_string(&mut contents).expect("Unable to read the Glimpse file");
  let glimpse: Glimpse = serde_json::from_str(&contents)?;
  Ok(glimpse)
}

#[derive(Serialize, Deserialize)]
pub struct NotGlimpse {
  #[serde(rename = "0")]
  d0: Vec<u64>,
  #[serde(rename = "1")]
  d1: Vec<u64>,
  #[serde(rename = "2")]
  d2: Vec<u64>,
  #[serde(rename = "3")]
  d3: Vec<u64>,
  #[serde(rename = "4")]
  d4: Vec<u64>,
  #[serde(rename = "5")]
  d5: Vec<u64>,
  #[serde(rename = "6")]
  d6: Vec<u64>,
  #[serde(rename = "7")]
  d7: Vec<u64>,
  #[serde(rename = "8")]
  d8: Vec<u64>,
  #[serde(rename = "9")]
  d9: Vec<u64>,
}

impl NotGlimpse {

  pub fn to_bmoc(&self) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(9,  self.d0.len() + self.d1.len() + self.d2.len() 
      + self.d3.len() + self.d4.len() + self.d5.len() + self.d6.len() + self.d7.len() 
      + self.d8.len() + self.d9.len());
    for h in &self.d0 {
      builder.push(0, *h, true);
    }
    for h in &self.d1 {
      builder.push(1, *h, true);
    }
    for h in &self.d2 {
      builder.push(2, *h, true);
    }
    for h in &self.d3 {
      builder.push(3, *h, true);
    }
    for h in &self.d4 {
      builder.push(4, *h, true);
    }
    for h in &self.d5 {
      builder.push(5, *h, true);
    }
    for h in &self.d6 {
      builder.push(6, *h, true);
    }
    for h in &self.d7 {
      builder.push(7, *h, true);
    }
    for h in &self.d8 {
      builder.push(8, *h, true);
    }
    for h in &self.d9 {
      builder.push(9, *h, true);
    }
    builder.to_bmoc_from_unordered()
  }
}

pub fn load_not_glimpse() -> Result<NotGlimpse> {
  let mut file = File::open("resources/glimpse_not.moc.json").expect("Unable to open the Not Glimpse file");
  let mut contents = String::new();
  file.read_to_string(&mut contents).expect("Unable to read the Glimpse file");
  let glimpse: NotGlimpse = serde_json::from_str(&contents)?;
  Ok(glimpse)
}




#[derive(Serialize, Deserialize)]
pub struct SandG {
  #[serde(rename = "5")]
  d5: Vec<u64>,
  #[serde(rename = "6")]
  d6: Vec<u64>,
  #[serde(rename = "7")]
  d7: Vec<u64>,
  #[serde(rename = "8")]
  d8: Vec<u64>,
  #[serde(rename = "9")]
  d9: Vec<u64>,
  #[serde(rename = "10")]
  d10: Vec<u64>,
  #[serde(rename = "11")]
  d11: Vec<u64>,
  #[serde(rename = "12")]
  d12: Vec<u64>,
}

impl SandG {
  pub fn to_bmoc(&self) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(12, self.d5.len() + self.d6.len() + self.d7.len()
      + self.d8.len() + self.d9.len() + self.d10.len() + self.d11.len() + self.d12.len());
    for h in &self.d5 {
      builder.push(5, *h, true);
    }
    for h in &self.d6 {
      builder.push(6, *h, true);
    }
    for h in &self.d7 {
      builder.push(7, *h, true);
    }
    for h in &self.d8 {
      builder.push(8, *h, true);
    }
    for h in &self.d9 {
      builder.push(9, *h, true);
    }
    for h in &self.d10 {
      builder.push(10, *h, true);
    }
    for h in &self.d11 {
      builder.push(11, *h, true);
    }
    for h in &self.d12 {
      builder.push(12, *h, true);
    }
    builder.to_bmoc_from_unordered()
  }
}

pub fn load_s_and_g() -> Result<SandG> {
  let mut file = File::open("resources/sdss_and_glimpse.moc.json").expect("Unable to open the s and g file");
  let mut contents = String::new();
  file.read_to_string(&mut contents).expect("Unable to read the s and g file");
  let s_and_g: SandG = serde_json::from_str(&contents)?;
  Ok(s_and_g)
}










#[cfg(test)]
mod tests {
  use super::*;
  
  #[test]
  fn test_sdss_not() {
    let sdss = load_sdss().unwrap().to_bmoc();
    
    /*println!("n_cell in 1: {}", sdss.d1.len());
    println!("n_cell in 2: {}", sdss.d2.len());
    println!("n_cell in 3: {}", sdss.d3.len());
    println!("n_cell in 4: {}", sdss.d4.len());
    println!("n_cell in 5: {}", sdss.d5.len());
    println!("n_cell in 6: {}", sdss.d6.len());
    println!("n_cell in 7: {}", sdss.d7.len());
    println!("n_cell in 8: {}", sdss.d8.len());
    println!("n_cell in 9: {}", sdss.d9.len());
    println!("n_cell in 10: {}", sdss.d10.len());
    println!("n_cell in 11: {}", sdss.d11.len());
    println!("n_cell in 12: {}", sdss.d12.len());
    let sdss_moc = sdss.to_bmoc();*/

    let sdss_not = load_sdss_not().unwrap().to_bmoc();
    
    let now = Instant::now();
    let not_sdss = sdss.not();
    println!("SDSS 'not' operation done in {} ms. In / Out sizes: {} / {}", 
             now.elapsed().as_millis(), sdss.entries.len(), not_sdss.entries.len());
    
    /*let l1 = sdss_not.entries.len();
    let l2 = not_sdss.entries.len();
    println!("l1: {}, l2: {}", l1, l2);
    for i in 0..l1.min(l2) {
      let c1 = sdss_not.from_raw_value(sdss_not.entries[i]);
      let c2 = not_sdss.from_raw_value(not_sdss.entries[i]);
      if c1.depth != c2.depth || c1.hash != c2.hash {
        println!("l: {:?}, r: {:?}, ln-1: {:?}, ln+1: {:?}, rn-1: {:?}, rn+1: {:?}", &c1, &c2,
                 sdss_not.from_raw_value(sdss_not.entries[i - 1]),
                 sdss_not.from_raw_value(sdss_not.entries[i + 1]),
                 not_sdss.from_raw_value(not_sdss.entries[i - 1]),
                 not_sdss.from_raw_value(not_sdss.entries[i + 1]));
        panic!("Toto");
      }
    }*/
    
    // sdss_not.assert_equals(&not_sdss);
    assert!(sdss_not.equals(&not_sdss));
  }

  #[test]
  fn test_glimpse_not() {
    let glimpse = load_glimpse().unwrap().to_bmoc();
    /*println!("n_cell in 4: {}", glimpse.d4.len());
    println!("n_cell in 5: {}", glimpse.d5.len());
    println!("n_cell in 6: {}", glimpse.d6.len());
    println!("n_cell in 7: {}", glimpse.d7.len());
    println!("n_cell in 8: {}", glimpse.d8.len());
    println!("n_cell in 9: {}", glimpse.d9.len());*/
    let not_glimpse = load_not_glimpse().unwrap().to_bmoc();
    
    let now = Instant::now();
    let glimpse_not = glimpse.not();
    println!("GLIMPSE 'not' operation done in {} ms. In / Out sizes: {} / {}", 
             now.elapsed().as_millis(), glimpse.entries.len(), not_glimpse.entries.len());
    
    assert!(not_glimpse.equals(&glimpse_not));
  }

  #[test]
  fn test_and() {
    let glimpse = load_glimpse().unwrap().to_bmoc();
    let sdss = load_sdss().unwrap().to_bmoc();

    let s_and_g = load_s_and_g().unwrap().to_bmoc();
    
    let now = Instant::now();
    let sandg = sdss.and(&glimpse);
    println!("SDSS/GLIMPSE 'and' operation done in {} ms. In / Out sizes: {} and {} => {}",
             now.elapsed().as_millis(), sdss.entries.len(), glimpse.entries.len(), sandg.entries.len());
    
    assert!(s_and_g.equals(&sandg));
  }

  
  #[test]
  fn test_or() {
    let glimpse = load_glimpse().unwrap().to_bmoc();
    let sdss = load_sdss().unwrap().to_bmoc();
    let s_or_g = load_s_or_g().unwrap().to_bmoc();
    
    let now = Instant::now();
    let sorg = &sdss.or(&glimpse);
    println!("SDSS/GLIMPSE 'or' operation done in {} ms. In / Out sizes: {} or {} => {}",
             now.elapsed().as_millis(), sdss.entries.len(), glimpse.entries.len(), sorg.entries.len());
    
    s_or_g.assert_equals(&sorg);
    assert!(s_or_g.equals(&sorg));
  }

  #[test]
  fn test_xor() {
    let glimpse = load_glimpse().unwrap().to_bmoc();
    let sdss = load_sdss().unwrap().to_bmoc();
    let s_xor_g = load_s_xor_g().unwrap().to_bmoc();
    
    let now = Instant::now();
    let sxorg = sdss.xor(&glimpse);
    println!("SDSS/GLIMPSE 'xor' operation done in {} ms. In / Out sizes: {} xor {} => {}",
             now.elapsed().as_millis(), sdss.entries.len(), glimpse.entries.len(), sxorg.entries.len());
    
    s_xor_g.assert_equals(&sxorg);
    assert!(s_xor_g.equals(&sdss.xor(&glimpse)));
  }

  #[test]
  fn test_sdss_build_from_ordered_input() {
    let sdss = load_sdss().unwrap().to_bmoc();
    let deep_size = sdss.deep_size();
    
    let mut builder = BMOCBuilderFixedDepth::new(sdss.get_depth_max(), true);

    let now = Instant::now();
    for h in sdss.flat_iter() {
      builder.push(h);
    }
    let sdss2 = builder.to_bmoc().unwrap();
    println!("SDSS 'build' operation (from SDSS flat iterator, deep size: {}) done in {} ms", deep_size, now.elapsed().as_millis());
    
    
    sdss.assert_equals(&sdss2);
    assert!(sdss.equals(&sdss2));
  }
  
}
