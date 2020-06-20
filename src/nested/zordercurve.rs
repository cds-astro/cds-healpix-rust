//! 2D Morton code implementation.
//! For BMI 2.0 references, see:
//! - [here](https://github.com/gnzlbg/bitwise/blob/master/src/word/morton.rs) to use BMI 2.0 pdep/pext instructions.
//! - [here](https://stackoverflow.com/questions/4909263/how-to-efficiently-de-interleave-bits-inverse-morton) another example with BMI2.0 instructions 
//! - [here](https://docs.rs/bitintr/0.2.0/src/bitintr/pdep.rs.html#68-87) code using BMI 2.0 and Rust macro 
//! - [The rust doc](https://doc.rust-lang.org/beta/core/arch/x86_64/) for the X86_64 architecture
//! - See also the comment of Julien Bilalte [here](https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/):
//! > If you can afford the luxury of compiling for BMI2-enabled CPUs, it just boils down to one 
//! > pdep / pext instruction per dimension`
//! > These instructions have latency/troughput of 3/1 on Haswell, so it ends up being much faster 
//! > than the LUT methods as well.
//! > (and also faster than a SIMD implementation that processes 4 morton codes at a time, 
//! > although one could probably avoid pipeline starvations by intertwining a SIMD implementation 
//! > with a few iterations of the pext/pdep versions to get even further speed out of it,
//! > if you’ve got the need to generate a stream of morton codes for some reason… :) )
// If needed one day, to encore 3D data (from one of the above ref, I think), from Julien Bilalte:
// encode3(u32 x, u32 y, u32 z) return _pdep_u32(z, 0x24924924) | _pdep_u32(y, 0x12492492) | _pdep_u32(x, 0x09249249);
// decode3(u32 code, u32 &outX, u32 &outY, u32 &outZ)
//  outX = _pext_u32(code, 0x09249249);
//  outY = _pext_u32(code, 0x12492492);
//  outZ = _pext_u32(code, 0x24924924);
use std::mem;

static LUPT_TO_HASH: [u16; 256] = [
  0x0000, 0x0001, 0x0004, 0x0005, 0x0010, 0x0011, 0x0014, 0x0015, 0x0040, 0x0041, 0x0044,
  0x0045, 0x0050, 0x0051, 0x0054, 0x0055, 0x0100, 0x0101, 0x0104, 0x0105, 0x0110, 0x0111,
  0x0114, 0x0115, 0x0140, 0x0141, 0x0144, 0x0145, 0x0150, 0x0151, 0x0154, 0x0155, 0x0400,
  0x0401, 0x0404, 0x0405, 0x0410, 0x0411, 0x0414, 0x0415, 0x0440, 0x0441, 0x0444, 0x0445,
  0x0450, 0x0451, 0x0454, 0x0455, 0x0500, 0x0501, 0x0504, 0x0505, 0x0510, 0x0511, 0x0514,
  0x0515, 0x0540, 0x0541, 0x0544, 0x0545, 0x0550, 0x0551, 0x0554, 0x0555, 0x1000, 0x1001,
  0x1004, 0x1005, 0x1010, 0x1011, 0x1014, 0x1015, 0x1040, 0x1041, 0x1044, 0x1045, 0x1050,
  0x1051, 0x1054, 0x1055, 0x1100, 0x1101, 0x1104, 0x1105, 0x1110, 0x1111, 0x1114, 0x1115,
  0x1140, 0x1141, 0x1144, 0x1145, 0x1150, 0x1151, 0x1154, 0x1155, 0x1400, 0x1401, 0x1404,
  0x1405, 0x1410, 0x1411, 0x1414, 0x1415, 0x1440, 0x1441, 0x1444, 0x1445, 0x1450, 0x1451,
  0x1454, 0x1455, 0x1500, 0x1501, 0x1504, 0x1505, 0x1510, 0x1511, 0x1514, 0x1515, 0x1540,
  0x1541, 0x1544, 0x1545, 0x1550, 0x1551, 0x1554, 0x1555, 0x4000, 0x4001, 0x4004, 0x4005,
  0x4010, 0x4011, 0x4014, 0x4015, 0x4040, 0x4041, 0x4044, 0x4045, 0x4050, 0x4051, 0x4054,
  0x4055, 0x4100, 0x4101, 0x4104, 0x4105, 0x4110, 0x4111, 0x4114, 0x4115, 0x4140, 0x4141,
  0x4144, 0x4145, 0x4150, 0x4151, 0x4154, 0x4155, 0x4400, 0x4401, 0x4404, 0x4405, 0x4410,
  0x4411, 0x4414, 0x4415, 0x4440, 0x4441, 0x4444, 0x4445, 0x4450, 0x4451, 0x4454, 0x4455,
  0x4500, 0x4501, 0x4504, 0x4505, 0x4510, 0x4511, 0x4514, 0x4515, 0x4540, 0x4541, 0x4544,
  0x4545, 0x4550, 0x4551, 0x4554, 0x4555, 0x5000, 0x5001, 0x5004, 0x5005, 0x5010, 0x5011,
  0x5014, 0x5015, 0x5040, 0x5041, 0x5044, 0x5045, 0x5050, 0x5051, 0x5054, 0x5055, 0x5100,
  0x5101, 0x5104, 0x5105, 0x5110, 0x5111, 0x5114, 0x5115, 0x5140, 0x5141, 0x5144, 0x5145,
  0x5150, 0x5151, 0x5154, 0x5155, 0x5400, 0x5401, 0x5404, 0x5405, 0x5410, 0x5411, 0x5414,
  0x5415, 0x5440, 0x5441, 0x5444, 0x5445, 0x5450, 0x5451, 0x5454, 0x5455, 0x5500, 0x5501,
  0x5504, 0x5505, 0x5510, 0x5511, 0x5514, 0x5515, 0x5540, 0x5541, 0x5544, 0x5545, 0x5550,
  0x5551, 0x5554, 0x5555];

static LUPT_TO_IJ_BYTE: [u16; 256] = [
  0x000, 0x001, 0x100, 0x101, 0x002, 0x003, 0x102, 0x103, 0x200, 0x201, 0x300, 0x301, 0x202,
  0x203, 0x302, 0x303, 0x004, 0x005, 0x104, 0x105, 0x006, 0x007, 0x106, 0x107, 0x204, 0x205,
  0x304, 0x305, 0x206, 0x207, 0x306, 0x307, 0x400, 0x401, 0x500, 0x501, 0x402, 0x403, 0x502,
  0x503, 0x600, 0x601, 0x700, 0x701, 0x602, 0x603, 0x702, 0x703, 0x404, 0x405, 0x504, 0x505,
  0x406, 0x407, 0x506, 0x507, 0x604, 0x605, 0x704, 0x705, 0x606, 0x607, 0x706, 0x707, 0x008,
  0x009, 0x108, 0x109, 0x00A, 0x00B, 0x10A, 0x10B, 0x208, 0x209, 0x308, 0x309, 0x20A, 0x20B,
  0x30A, 0x30B, 0x00C, 0x00D, 0x10C, 0x10D, 0x00E, 0x00F, 0x10E, 0x10F, 0x20C, 0x20D, 0x30C,
  0x30D, 0x20E, 0x20F, 0x30E, 0x30F, 0x408, 0x409, 0x508, 0x509, 0x40A, 0x40B, 0x50A, 0x50B,
  0x608, 0x609, 0x708, 0x709, 0x60A, 0x60B, 0x70A, 0x70B, 0x40C, 0x40D, 0x50C, 0x50D, 0x40E,
  0x40F, 0x50E, 0x50F, 0x60C, 0x60D, 0x70C, 0x70D, 0x60E, 0x60F, 0x70E, 0x70F, 0x800, 0x801,
  0x900, 0x901, 0x802, 0x803, 0x902, 0x903, 0xA00, 0xA01, 0xB00, 0xB01, 0xA02, 0xA03, 0xB02,
  0xB03, 0x804, 0x805, 0x904, 0x905, 0x806, 0x807, 0x906, 0x907, 0xA04, 0xA05, 0xB04, 0xB05,
  0xA06, 0xA07, 0xB06, 0xB07, 0xC00, 0xC01, 0xD00, 0xD01, 0xC02, 0xC03, 0xD02, 0xD03, 0xE00,
  0xE01, 0xF00, 0xF01, 0xE02, 0xE03, 0xF02, 0xF03, 0xC04, 0xC05, 0xD04, 0xD05, 0xC06, 0xC07,
  0xD06, 0xD07, 0xE04, 0xE05, 0xF04, 0xF05, 0xE06, 0xE07, 0xF06, 0xF07, 0x808, 0x809, 0x908,
  0x909, 0x80A, 0x80B, 0x90A, 0x90B, 0xA08, 0xA09, 0xB08, 0xB09, 0xA0A, 0xA0B, 0xB0A, 0xB0B,
  0x80C, 0x80D, 0x90C, 0x90D, 0x80E, 0x80F, 0x90E, 0x90F, 0xA0C, 0xA0D, 0xB0C, 0xB0D, 0xA0E,
  0xA0F, 0xB0E, 0xB0F, 0xC08, 0xC09, 0xD08, 0xD09, 0xC0A, 0xC0B, 0xD0A, 0xD0B, 0xE08, 0xE09,
  0xF08, 0xF09, 0xE0A, 0xE0B, 0xF0A, 0xF0B, 0xC0C, 0xC0D, 0xD0C, 0xD0D, 0xC0E, 0xC0F, 0xD0E,
  0xD0F, 0xE0C, 0xE0D, 0xF0C, 0xF0D, 0xE0E, 0xE0F, 0xF0E, 0xF0F];

static LUPT_TO_IJ_SHORT: [u32; 256] = [
  0x00000, 0x00001, 0x10000, 0x10001, 0x00002, 0x00003, 0x10002, 0x10003, 0x20000, 0x20001,
  0x30000, 0x30001, 0x20002, 0x20003, 0x30002, 0x30003, 0x00004, 0x00005, 0x10004, 0x10005,
  0x00006, 0x00007, 0x10006, 0x10007, 0x20004, 0x20005, 0x30004, 0x30005, 0x20006, 0x20007,
  0x30006, 0x30007, 0x40000, 0x40001, 0x50000, 0x50001, 0x40002, 0x40003, 0x50002, 0x50003,
  0x60000, 0x60001, 0x70000, 0x70001, 0x60002, 0x60003, 0x70002, 0x70003, 0x40004, 0x40005,
  0x50004, 0x50005, 0x40006, 0x40007, 0x50006, 0x50007, 0x60004, 0x60005, 0x70004, 0x70005,
  0x60006, 0x60007, 0x70006, 0x70007, 0x00008, 0x00009, 0x10008, 0x10009, 0x0000A, 0x0000B,
  0x1000A, 0x1000B, 0x20008, 0x20009, 0x30008, 0x30009, 0x2000A, 0x2000B, 0x3000A, 0x3000B,
  0x0000C, 0x0000D, 0x1000C, 0x1000D, 0x0000E, 0x0000F, 0x1000E, 0x1000F, 0x2000C, 0x2000D,
  0x3000C, 0x3000D, 0x2000E, 0x2000F, 0x3000E, 0x3000F, 0x40008, 0x40009, 0x50008, 0x50009,
  0x4000A, 0x4000B, 0x5000A, 0x5000B, 0x60008, 0x60009, 0x70008, 0x70009, 0x6000A, 0x6000B,
  0x7000A, 0x7000B, 0x4000C, 0x4000D, 0x5000C, 0x5000D, 0x4000E, 0x4000F, 0x5000E, 0x5000F,
  0x6000C, 0x6000D, 0x7000C, 0x7000D, 0x6000E, 0x6000F, 0x7000E, 0x7000F, 0x80000, 0x80001,
  0x90000, 0x90001, 0x80002, 0x80003, 0x90002, 0x90003, 0xA0000, 0xA0001, 0xB0000, 0xB0001,
  0xA0002, 0xA0003, 0xB0002, 0xB0003, 0x80004, 0x80005, 0x90004, 0x90005, 0x80006, 0x80007,
  0x90006, 0x90007, 0xA0004, 0xA0005, 0xB0004, 0xB0005, 0xA0006, 0xA0007, 0xB0006, 0xB0007,
  0xC0000, 0xC0001, 0xD0000, 0xD0001, 0xC0002, 0xC0003, 0xD0002, 0xD0003, 0xE0000, 0xE0001,
  0xF0000, 0xF0001, 0xE0002, 0xE0003, 0xF0002, 0xF0003, 0xC0004, 0xC0005, 0xD0004, 0xD0005,
  0xC0006, 0xC0007, 0xD0006, 0xD0007, 0xE0004, 0xE0005, 0xF0004, 0xF0005, 0xE0006, 0xE0007,
  0xF0006, 0xF0007, 0x80008, 0x80009, 0x90008, 0x90009, 0x8000A, 0x8000B, 0x9000A, 0x9000B,
  0xA0008, 0xA0009, 0xB0008, 0xB0009, 0xA000A, 0xA000B, 0xB000A, 0xB000B, 0x8000C, 0x8000D,
  0x9000C, 0x9000D, 0x8000E, 0x8000F, 0x9000E, 0x9000F, 0xA000C, 0xA000D, 0xB000C, 0xB000D,
  0xA000E, 0xA000F, 0xB000E, 0xB000F, 0xC0008, 0xC0009, 0xD0008, 0xD0009, 0xC000A, 0xC000B,
  0xD000A, 0xD000B, 0xE0008, 0xE0009, 0xF0008, 0xF0009, 0xE000A, 0xE000B, 0xF000A, 0xF000B,
  0xC000C, 0xC000D, 0xD000C, 0xD000D, 0xC000E, 0xC000F, 0xD000E, 0xD000F, 0xE000C, 0xE000D,
  0xF000C, 0xF000D, 0xE000E, 0xE000F, 0xF000E, 0xF000F];

static LUPT_TO_IJ_INT: [u64; 256] = [
  0x000000000, 0x000000001, 0x100000000, 0x100000001, 0x000000002, 0x000000003,
  0x100000002, 0x100000003, 0x200000000, 0x200000001, 0x300000000, 0x300000001,
  0x200000002, 0x200000003, 0x300000002, 0x300000003, 0x000000004, 0x000000005,
  0x100000004, 0x100000005, 0x000000006, 0x000000007, 0x100000006, 0x100000007,
  0x200000004, 0x200000005, 0x300000004, 0x300000005, 0x200000006, 0x200000007,
  0x300000006, 0x300000007, 0x400000000, 0x400000001, 0x500000000, 0x500000001,
  0x400000002, 0x400000003, 0x500000002, 0x500000003, 0x600000000, 0x600000001,
  0x700000000, 0x700000001, 0x600000002, 0x600000003, 0x700000002, 0x700000003,
  0x400000004, 0x400000005, 0x500000004, 0x500000005, 0x400000006, 0x400000007,
  0x500000006, 0x500000007, 0x600000004, 0x600000005, 0x700000004, 0x700000005,
  0x600000006, 0x600000007, 0x700000006, 0x700000007, 0x000000008, 0x000000009,
  0x100000008, 0x100000009, 0x00000000A, 0x00000000B, 0x10000000A, 0x10000000B,
  0x200000008, 0x200000009, 0x300000008, 0x300000009, 0x20000000A, 0x20000000B,
  0x30000000A, 0x30000000B, 0x00000000C, 0x00000000D, 0x10000000C, 0x10000000D,
  0x00000000E, 0x00000000F, 0x10000000E, 0x10000000F, 0x20000000C, 0x20000000D,
  0x30000000C, 0x30000000D, 0x20000000E, 0x20000000F, 0x30000000E, 0x30000000F,
  0x400000008, 0x400000009, 0x500000008, 0x500000009, 0x40000000A, 0x40000000B,
  0x50000000A, 0x50000000B, 0x600000008, 0x600000009, 0x700000008, 0x700000009,
  0x60000000A, 0x60000000B, 0x70000000A, 0x70000000B, 0x40000000C, 0x40000000D,
  0x50000000C, 0x50000000D, 0x40000000E, 0x40000000F, 0x50000000E, 0x50000000F,
  0x60000000C, 0x60000000D, 0x70000000C, 0x70000000D, 0x60000000E, 0x60000000F,
  0x70000000E, 0x70000000F, 0x800000000, 0x800000001, 0x900000000, 0x900000001,
  0x800000002, 0x800000003, 0x900000002, 0x900000003, 0xA00000000, 0xA00000001,
  0xB00000000, 0xB00000001, 0xA00000002, 0xA00000003, 0xB00000002, 0xB00000003,
  0x800000004, 0x800000005, 0x900000004, 0x900000005, 0x800000006, 0x800000007,
  0x900000006, 0x900000007, 0xA00000004, 0xA00000005, 0xB00000004, 0xB00000005,
  0xA00000006, 0xA00000007, 0xB00000006, 0xB00000007, 0xC00000000, 0xC00000001,
  0xD00000000, 0xD00000001, 0xC00000002, 0xC00000003, 0xD00000002, 0xD00000003,
  0xE00000000, 0xE00000001, 0xF00000000, 0xF00000001, 0xE00000002, 0xE00000003,
  0xF00000002, 0xF00000003, 0xC00000004, 0xC00000005, 0xD00000004, 0xD00000005,
  0xC00000006, 0xC00000007, 0xD00000006, 0xD00000007, 0xE00000004, 0xE00000005,
  0xF00000004, 0xF00000005, 0xE00000006, 0xE00000007, 0xF00000006, 0xF00000007,
  0x800000008, 0x800000009, 0x900000008, 0x900000009, 0x80000000A, 0x80000000B,
  0x90000000A, 0x90000000B, 0xA00000008, 0xA00000009, 0xB00000008, 0xB00000009,
  0xA0000000A, 0xA0000000B, 0xB0000000A, 0xB0000000B, 0x80000000C, 0x80000000D,
  0x90000000C, 0x90000000D, 0x80000000E, 0x80000000F, 0x90000000E, 0x90000000F,
  0xA0000000C, 0xA0000000D, 0xB0000000C, 0xB0000000D, 0xA0000000E, 0xA0000000F,
  0xB0000000E, 0xB0000000F, 0xC00000008, 0xC00000009, 0xD00000008, 0xD00000009,
  0xC0000000A, 0xC0000000B, 0xD0000000A, 0xD0000000B, 0xE00000008, 0xE00000009,
  0xF00000008, 0xF00000009, 0xE0000000A, 0xE0000000B, 0xF0000000A, 0xF0000000B,
  0xC0000000C, 0xC0000000D, 0xD0000000C, 0xD0000000D, 0xC0000000E, 0xC0000000F,
  0xD0000000E, 0xD0000000F, 0xE0000000C, 0xE0000000D, 0xF0000000C, 0xF0000000D,
  0xE0000000E, 0xE0000000F, 0xF0000000E, 0xF0000000F];

pub trait ZOrderCurve: Sync + Send {

  fn ij2h(&self, i: u32, j: u32) -> u64 {
    self.i02h(i) | self.oj2h(j)
  }

  fn i02h(&self, i: u32) -> u64;
  fn oj2h(&self, j: u32) -> u64 {
    self.i02h(j) << 1
  }
  fn h2ij(&self, h: u64) -> u64;

  fn h2i0(&self, h: u64) -> u64 {
    self.h2ij(h)
  }

  fn ij2i(&self, ij: u64) -> u32;
  fn ij2j(&self, ij: u64) -> u32;

  fn xy2h(&self, x: f64, y: f64) -> u64 {
    self.ij2h(x as u32, y as u32)
  }

}

struct EmptyZOC;
impl ZOrderCurve for EmptyZOC {
  fn ij2h(&self, _i: u32, _j: u32) -> u64 { 0 }
  fn i02h(&self, _i: u32) -> u64 { 0 }
  fn h2ij(&self, _h: u64) -> u64 { 0 }
  fn h2i0(&self, _h: u64) -> u64 { 0 }
  fn ij2i(&self, _ij: u64) -> u32 { 0 }
  fn ij2j(&self, _ij: u64) -> u32 { 0 }
}

struct SmallZOC;
impl ZOrderCurve for SmallZOC {
  fn i02h(&self, i: u32) -> u64 {
    LUPT_TO_HASH[(i as u8) as usize] as u64
  }
  fn h2ij(&self, h: u64) -> u64 {
    let bytes:  [u8; 2] = (h as u16).to_le_bytes();
    (   LUPT_TO_IJ_BYTE[bytes[0] as usize]
      | LUPT_TO_IJ_BYTE[bytes[1] as usize] << 4
    ) as u64
  }
  fn ij2i(&self, ij: u64) -> u32 {
    (ij as u32) & 0x000000FF
  }
  fn ij2j(&self, ij: u64) -> u32 {
    (ij as u32) >> 8
  }
}

struct MediuZOC;
impl ZOrderCurve for MediuZOC {
  fn i02h(&self, i: u32) -> u64 {
    let bytes: [u8; 2] = (i as u16).to_le_bytes();
       LUPT_TO_HASH[bytes[0] as usize] as u64
    | (LUPT_TO_HASH[bytes[1] as usize] as u64) << 16
  }
  fn h2ij(&self, h: u64) -> u64 {
    let bytes: [u8; 4] = (h as u32).to_le_bytes();
    (   LUPT_TO_IJ_SHORT[bytes[0] as usize]
      | LUPT_TO_IJ_SHORT[bytes[1] as usize] <<  4
      | LUPT_TO_IJ_SHORT[bytes[2] as usize] <<  8
      | LUPT_TO_IJ_SHORT[bytes[3] as usize] << 12
    ) as u64
  }
  fn ij2i(&self, ij: u64) -> u32 {
    (ij as u32) & 0x0000FFFF
  }
  fn ij2j(&self, ij: u64) -> u32 {
    (ij as u32) >> 16
  }
}

pub struct LargeZOC;
impl ZOrderCurve for LargeZOC {
  fn i02h(&self, i: u32) -> u64 {
    // to/from_le do nothing on x86 architectures (which are in LE), so not perf penalty
    let bytes: [u8; 4] = i.to_le_bytes();
    /*(     LUPT_TO_HASH[bytes[0] as usize] as u64
				| (LUPT_TO_HASH[bytes[1] as usize] as u64) << 16
				| (LUPT_TO_HASH[bytes[2] as usize] as u64) << 32
				| (LUPT_TO_HASH[bytes[3] as usize] as u64) << 48
		) as u64*/
    // Portability on BE architecture to be tested!! (if not portable, use the above formulation)
    u64::from_le(
      unsafe {
        mem::transmute::<[u16; 4], u64>([
          LUPT_TO_HASH[bytes[0] as usize],
          LUPT_TO_HASH[bytes[1] as usize],
          LUPT_TO_HASH[bytes[2] as usize],
          LUPT_TO_HASH[bytes[3] as usize]])
      }
    )
  }
  fn h2ij(&self, h: u64) -> u64 {
    let bytes: [u8; 8] = h.to_le_bytes();
      LUPT_TO_IJ_INT[bytes[0] as usize] as u64
    | LUPT_TO_IJ_INT[bytes[1] as usize] <<  4
    | LUPT_TO_IJ_INT[bytes[2] as usize] <<  8
    | LUPT_TO_IJ_INT[bytes[3] as usize] << 12
    | LUPT_TO_IJ_INT[bytes[4] as usize] << 16
    | LUPT_TO_IJ_INT[bytes[5] as usize] << 20
    | LUPT_TO_IJ_INT[bytes[6] as usize] << 24
    | LUPT_TO_IJ_INT[bytes[7] as usize] << 28
  }
  fn ij2i(&self, ij: u64) -> u32 {
    ij as u32
  }
  fn ij2j(&self, ij: u64) -> u32 {
    (ij >> 32) as u32
  }
}

#[cfg(test)]
struct SmallZOCxor;
#[cfg(test)]
impl ZOrderCurve for SmallZOCxor {
  fn ij2h(&self, i: u32, j: u32) -> u64 {
    let mut i = i as u16;
    let mut j = j as u16;
    i |= j << 8;
    j = (i ^ (i >> 4)) & 0x00F0u16; i = i ^ j ^ (j << 4);
    j = (i ^ (i >> 2)) & 0x0C0Cu16; i = i ^ j ^ (j << 2);
    j = (i ^ (i >> 1)) & 0x2222u16; i = i ^ j ^ (j << 1);
    i as u64
  }
  fn i02h(&self, i: u32) -> u64 {
    let mut i = i as u16;
    i = ((i << 4) | i) & 0x0F0Fu16;
    i = ((i << 2) | i) & 0x3333u16;
    i = ((i << 1) | i) & 0x5555u16;
    i as u64
  }
  fn h2ij(&self, h: u64) -> u64 {
    let mut h = h as u16;
    let mut t
      = (h ^ (h >> 1)) & 0x2222u16; h = h ^ t ^ (t << 1);
    t = (h ^ (h >> 2)) & 0x0C0Cu16; h = h ^ t ^ (t << 2);
    t = (h ^ (h >> 4)) & 0x00F0u16; h = h ^ t ^ (t << 4);
    h as u64
  }
  fn h2i0(&self, h: u64) -> u64 {
    let mut h = h as u16;
    h = ((h >> 1) | h) & 0x3333u16;
    h = ((h >> 2) | h) & 0x0F0Fu16;
    h = ((h >> 4) | h) & 0x00FFu16;
    h as u64
  }
  fn ij2i(&self, ij: u64) -> u32 {
    (ij as u32) & 0x000000FF
  }
  fn ij2j(&self, ij: u64) -> u32 {
    (ij as u32) >> 8
  }
}

#[cfg(test)]
struct MediuZOCxor;
#[cfg(test)]
impl ZOrderCurve for MediuZOCxor {
  fn ij2h(&self, mut i: u32, mut j: u32) -> u64 {
    i |= j << 16;
    j = (i ^ (i >> 8)) & 0x0000FF00u32; i = i ^ j ^ (j << 8);
    j = (i ^ (i >> 4)) & 0x00F000F0u32; i = i ^ j ^ (j << 4);
    j = (i ^ (i >> 2)) & 0x0C0C0C0Cu32; i = i ^ j ^ (j << 2);
    j = (i ^ (i >> 1)) & 0x22222222u32; i = i ^ j ^ (j << 1);
    i as u64
  }
  fn i02h(&self, mut i: u32) -> u64 {
    i = ((i << 8) | i) & 0x00FF00FFu32;
    i = ((i << 4) | i) & 0x0F0F0F0Fu32;
    i = ((i << 2) | i) & 0x33333333u32;
    i = ((i << 1) | i) & 0x55555555u32;
    i as u64
  }
  fn h2ij(&self, h: u64) -> u64 {
    let mut h = h as u32;
    let mut t
      = (h ^ (h >> 1)) & 0x22222222u32; h = h ^ t ^ (t << 1);
    t = (h ^ (h >> 2)) & 0x0C0C0C0Cu32; h = h ^ t ^ (t << 2);
    t = (h ^ (h >> 4)) & 0x00F000F0u32; h = h ^ t ^ (t << 4);
    t = (h ^ (h >> 8)) & 0x0000FF00u32; h = h ^ t ^ (t << 8);
    h as u64
  }
  fn h2i0(&self, mut h: u64) -> u64 {
    h = ((h >> 1) | h) & 0x33333333u64;
    h = ((h >> 2) | h) & 0x0F0F0F0Fu64;
    h = ((h >> 4) | h) & 0x00FF00FFu64;
    h = ((h >> 8) | h) & 0x0000FFFFu64;
    h
  }
  fn ij2i(&self, ij: u64) -> u32 {
    (ij as u32) & 0x0000FFFF
  }
  fn ij2j(&self, ij: u64) -> u32 { (ij as u32) >> 16 }
}

// #[cfg(any(test, bench))]
pub struct LargeZOCxor;
// #[cfg(any(test, bench))]
impl ZOrderCurve for LargeZOCxor {
  fn ij2h(&self, i: u32, j: u32) -> u64 {
    let mut h = ((j as u64) << 32) | (i as u64);
    let mut t
      = (h ^ (h >> 16)) & 0x00000000FFFF0000u64; h = h ^ t ^ (t << 16);
    t = (h ^ (h >>  8)) & 0x0000FF000000FF00u64; h = h ^ t ^ (t <<  8);
    t = (h ^ (h >>  4)) & 0x00F000F000F000F0u64; h = h ^ t ^ (t <<  4);
    t = (h ^ (h >>  2)) & 0x0C0C0C0C0C0C0C0Cu64; h = h ^ t ^ (t <<  2);
    t = (h ^ (h >>  1)) & 0x2222222222222222u64; h = h ^ t ^ (t <<  1);
    h
  }
  fn i02h(&self, i: u32) -> u64 {
    let mut i = i as u64;
    i = ((i << 16) | i) & 0x0000FFFF0000FFFFu64;
    i = ((i <<  8) | i) & 0x00FF00FF00FF00FFu64;
    i = ((i <<  4) | i) & 0x0F0F0F0F0F0F0F0Fu64;
    i = ((i <<  2) | i) & 0x3333333333333333u64;
    i = ((i <<  1) | i) & 0x5555555555555555u64;
    i
  }
  fn h2ij(&self, mut h: u64) -> u64 {
    let mut t
      = (h ^ (h >>  1)) & 0x2222222222222222u64; h = h ^ t ^ (t <<  1);
    t = (h ^ (h >>  2)) & 0x0C0C0C0C0C0C0C0Cu64; h = h ^ t ^ (t <<  2);
    t = (h ^ (h >>  4)) & 0x00F000F000F000F0u64; h = h ^ t ^ (t <<  4);
    t = (h ^ (h >>  8)) & 0x0000FF000000FF00u64; h = h ^ t ^ (t <<  8);
    t = (h ^ (h >> 16)) & 0x00000000FFFF0000u64; h = h ^ t ^ (t << 16);
    h
  }
  fn h2i0(&self, mut h: u64) -> u64 {
    h = ((h >>  1) | h) & 0x3333333333333333u64;
    h = ((h >>  2) | h) & 0x0F0F0F0F0F0F0F0Fu64;
    h = ((h >>  4) | h) & 0x00FF00FF00FF00FFu64;
    h = ((h >>  8) | h) & 0x0000FFFF0000FFFFu64;
    h = ((h >> 16) | h) & 0x00000000FFFFFFFFu64;
    h
  }
  fn ij2i(&self, ij: u64) -> u32 {
    ij as u32
  }
  fn ij2j(&self, ij: u64) -> u32 { (ij >> 32) as u32 }
}

// Instead of multiplying the #[cfg()], I should learn how to write macros... 

#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
pub struct SmallZOCbmi;
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
impl ZOrderCurve for SmallZOCbmi {
  
  fn ij2h(&self, i: u32, j: u32) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pdep_u32;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pdep_u32;
    unsafe {
      (_pdep_u32(i, 0x00005555u32) | _pdep_u32(j, 0x0000AAAAu32)) as u64
    }
  }
  fn i02h(&self, i: u32) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pdep_u32;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pdep_u32;
    unsafe {
      _pdep_u32(i, 0x00005555u32) as u64
    }
  }
  fn oj2h(&self, j: u32) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pdep_u32;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pdep_u32;
    unsafe {
      _pdep_u32(j, 0x0000AAAAu32) as u64
    }
  }
  fn h2ij(&self, h: u64) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pext_u32;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pext_u32;
    unsafe {
      (_pext_u32(h as u32, 0x00005555u32) | (_pext_u32(h as u32, 0x0000AAAAu32) << 8)) as u64
    }
  }
  fn h2i0(&self, h: u64) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pext_u32;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pext_u32;
    unsafe {
      _pext_u32(h as u32, 0x00005555u32) as u64
    }
  }
  fn ij2i(&self, ij: u64) -> u32 { (ij as u32) & 0x000000FFu32 }
  fn ij2j(&self, ij: u64) -> u32 {
    (ij as u32) >> 8
  }
}

#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
struct MediuZOCbmi;
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
impl ZOrderCurve for MediuZOCbmi {
  
  fn ij2h(&self, i: u32, j: u32) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pdep_u32;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pdep_u32;
    unsafe {
      (_pdep_u32(i, 0x55555555u32) | _pdep_u32(j, 0xAAAAAAAAu32)) as u64
    }
  }
  fn i02h(&self, i: u32) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pdep_u32;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pdep_u32;
    unsafe {
      _pdep_u32(i, 0x55555555u32) as u64
    }
  }
  fn oj2h(&self, i: u32) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pdep_u32;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pdep_u32;
    unsafe {
      _pdep_u32(i, 0xAAAAAAAAu32) as u64
    }
  }
  fn h2ij(&self, h: u64) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pext_u32;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pext_u32;
    unsafe {
      (_pext_u32(h as u32, 0x55555555u32) | (_pext_u32(h as u32, 0xAAAAAAAAu32) << 16)) as u64
    }
  }
  fn h2i0(&self, h: u64) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pext_u32;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pext_u32;
    unsafe {
      _pext_u32(h as u32, 0x55555555u32) as u64
    }
  }
  fn ij2i(&self, ij: u64) -> u32 { (ij as u32) & 0x0000FFFFu32 }
  fn ij2j(&self, ij: u64) -> u32 { (ij as u32) >> 16 }
}

#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
pub struct LargeZOCbmi;
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
impl ZOrderCurve for LargeZOCbmi {
  
  fn ij2h(&self, i: u32, j: u32) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pdep_u32;
    #[cfg(target_arch = "x86")]
    unsafe {
          (_pdep_u32(i & 0x0000FFFFu32, 0x55555555u32) as u64) 
        | (_pdep_u32(j & 0x0000FFFFu32, 0xAAAAAAAAu32) as u64) 
        | (_pdep_u32(i >> 16, 0x55555555u32) as u64) << 32
        | (_pdep_u32(j >> 16, 0xAAAAAAAAu32) as u64) << 32
    }
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pdep_u64;
    #[cfg(target_arch = "x86_64")]
    unsafe {
        _pdep_u64(i as u64, 0x5555555555555555u64) | _pdep_u64(j as u64, 0xAAAAAAAAAAAAAAAAu64)
    }
  }
  fn i02h(&self, i: u32) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pdep_u32;
    #[cfg(target_arch = "x86")]
    unsafe {
          (_pdep_u32(i & 0x0000FFFFu32, 0x55555555u32) as u64)
        | (_pdep_u32(i >> 16, 0x55555555u32) as u64) << 32
    }
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pdep_u64;
    #[cfg(target_arch = "x86_64")]
    unsafe {
      _pdep_u64(i as u64, 0x5555555555555555u64)
    }
  }
  fn oj2h(&self, i: u32) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pdep_u32;
    #[cfg(target_arch = "x86")]
    unsafe {
      (_pdep_u32(i & 0x0000FFFFu32, 0xAAAAAAAAu32) as u64)
        | (_pdep_u32(i >> 16, 0xAAAAAAAAu32) as u64) << 32
    }
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pdep_u64;
    #[cfg(target_arch = "x86_64")]
      unsafe {
        _pdep_u64(i as u64, 0xAAAAAAAAAAAAAAAAu64)
      }
  }
  fn h2ij(&self, h: u64) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pext_u32;
    #[cfg(target_arch = "x86")]
    unsafe {
          (_pext_u32(h as u32, 0x55555555u32) as u64)
        | (_pext_u32((h >> 32) as u32, 0x55555555u32) as u64) << 16
        | (_pext_u32(h as u32, 0xAAAAAAAAu32) as u64) << 32
        | (_pext_u32((h >> 32) as u32, 0xAAAAAAAAu32) as u64) << 48
    }
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pext_u64;
    #[cfg(target_arch = "x86_64")]
    unsafe {
      _pext_u64(h, 0x5555555555555555u64) | (_pext_u64(h, 0xAAAAAAAAAAAAAAAAu64) << 32)
    }
  }
  fn h2i0(&self, h: u64) -> u64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::_pext_u32;
    #[cfg(target_arch = "x86")]
    unsafe {
          (_pext_u32(h as u32, 0x55555555u32) as u64) 
        | (_pext_u32((h >> 32) as u32, 0x55555555u32) as u64) << 16
    }
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::_pext_u64;
    #[cfg(target_arch = "x86_64")]
    unsafe {
      _pext_u64(h, 0x5555555555555555u64)
    }
  }
  fn ij2i(&self, ij: u64) -> u32 {
    ij as u32
  }
  fn ij2j(&self, ij: u64) -> u32 { (ij >> 32) as u32 }
}


static EMPTY_ZOC: EmptyZOC = EmptyZOC;

static SMALL_ZOC_LUT: SmallZOC = SmallZOC;
static MEDIU_ZOC_LUT: MediuZOC = MediuZOC;
pub static LARGE_ZOC_LUT: LargeZOC = LargeZOC;

#[cfg(test)]
static SMALL_ZOC_XOR: SmallZOCxor = SmallZOCxor;
#[cfg(test)]
static MEDIU_ZOC_XOR: MediuZOCxor = MediuZOCxor;
// #[cfg(any(test, bench))]
pub static LARGE_ZOC_XOR: LargeZOCxor = LargeZOCxor;

#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
static SMALL_ZOC_BMI: SmallZOCbmi = SmallZOCbmi;
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
static MEDIU_ZOC_BMI: MediuZOCbmi = MediuZOCbmi;
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
pub static LARGE_ZOC_BMI: LargeZOCbmi = LargeZOCbmi;

pub enum ZOC {
  EMPTY,  // for depth = 0
  SMALL,  // for depth in [1, 8]
  MEDIUM, // for depth in [9, 16]
  LARGE,  // for depth in [17, 29]
}

#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
impl ZOrderCurve for ZOC {
  fn ij2h(&self, i: u32, j: u32) -> u64 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.ij2h(i, j),
      ZOC::SMALL => SMALL_ZOC_BMI.ij2h(i, j),
      ZOC::MEDIUM => MEDIU_ZOC_BMI.ij2h(i, j),
      ZOC::LARGE => LARGE_ZOC_BMI.ij2h(i, j),
    }
  }
  fn i02h(&self, i: u32) -> u64 {
    match self {
      ZOC:: EMPTY => EMPTY_ZOC.i02h(i),
      ZOC::SMALL => SMALL_ZOC_BMI.i02h(i),
      ZOC::MEDIUM => MEDIU_ZOC_BMI.i02h(i),
      ZOC::LARGE => LARGE_ZOC_BMI.i02h(i),
    }
  }
  fn oj2h(&self, j: u32) -> u64 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.oj2h(j),
      ZOC::SMALL => SMALL_ZOC_BMI.oj2h(j),
      ZOC::MEDIUM => MEDIU_ZOC_BMI.oj2h(j),
      ZOC::LARGE => LARGE_ZOC_BMI.oj2h(j),
    }
  }
  fn h2ij(&self, h: u64) -> u64 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.h2ij(h),
      ZOC::SMALL => SMALL_ZOC_BMI.h2ij(h),
      ZOC::MEDIUM => MEDIU_ZOC_BMI.h2ij(h),
      ZOC::LARGE => LARGE_ZOC_BMI.h2ij(h),
    }
  }
  fn h2i0(&self, h: u64) -> u64 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.h2i0(h),
      ZOC::SMALL => SMALL_ZOC_BMI.h2i0(h),
      ZOC::MEDIUM => MEDIU_ZOC_BMI.h2i0(h),
      ZOC::LARGE => LARGE_ZOC_BMI.h2i0(h),
    }
  }
  fn ij2i(&self, ij: u64) -> u32 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.ij2i(ij),
      ZOC::SMALL => SMALL_ZOC_BMI.ij2i(ij),
      ZOC::MEDIUM => MEDIU_ZOC_BMI.ij2i(ij),
      ZOC::LARGE => LARGE_ZOC_BMI.ij2i(ij),
    }
  }
  fn ij2j(&self, ij: u64) -> u32 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.ij2j(ij),
      ZOC::SMALL => SMALL_ZOC_BMI.ij2j(ij),
      ZOC::MEDIUM => MEDIU_ZOC_BMI.ij2j(ij),
      ZOC::LARGE => LARGE_ZOC_BMI.ij2j(ij),
    }
  }
  fn xy2h(&self, x: f64, y: f64) -> u64{
    match self {
      ZOC::EMPTY => EMPTY_ZOC.xy2h(x, y),
      ZOC::SMALL => SMALL_ZOC_BMI.xy2h(x, y),
      ZOC::MEDIUM => MEDIU_ZOC_BMI.xy2h(x, y),
      ZOC::LARGE => LARGE_ZOC_BMI.xy2h(x, y),
    }
  }
}

#[cfg(not(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2")))]
impl ZOrderCurve for ZOC {
  fn ij2h(&self, i: u32, j: u32) -> u64 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.ij2h(i, j),
      ZOC::SMALL => SMALL_ZOC_LUT.ij2h(i, j),
      ZOC::MEDIUM => MEDIU_ZOC_LUT.ij2h(i, j),
      ZOC::LARGE => LARGE_ZOC_LUT.ij2h(i, j),
    }
  }
  fn i02h(&self, i: u32) -> u64 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.i02h(i),
      ZOC::SMALL => SMALL_ZOC_LUT.i02h(i),
      ZOC::MEDIUM => MEDIU_ZOC_LUT.i02h(i),
      ZOC::LARGE => LARGE_ZOC_LUT.i02h(i),
    }
  }
  fn oj2h(&self, j: u32) -> u64 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.oj2h(j),
      ZOC::SMALL => SMALL_ZOC_LUT.oj2h(j),
      ZOC::MEDIUM => MEDIU_ZOC_LUT.oj2h(j),
      ZOC::LARGE => LARGE_ZOC_LUT.oj2h(j),
    }
  }
  fn h2ij(&self, h: u64) -> u64 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.h2ij(h),
      ZOC::SMALL => SMALL_ZOC_LUT.h2ij(h),
      ZOC::MEDIUM => MEDIU_ZOC_LUT.h2ij(h),
      ZOC::LARGE => LARGE_ZOC_LUT.h2ij(h),
    }
  }
  fn h2i0(&self, h: u64) -> u64 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.h2i0(h),
      ZOC::SMALL => SMALL_ZOC_LUT.h2i0(h),
      ZOC::MEDIUM => MEDIU_ZOC_LUT.h2i0(h),
      ZOC::LARGE => LARGE_ZOC_LUT.h2i0(h),
    }
  }
  fn ij2i(&self, ij: u64) -> u32 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.ij2i(ij),
      ZOC::SMALL => SMALL_ZOC_LUT.ij2i(ij),
      ZOC::MEDIUM => MEDIU_ZOC_LUT.ij2i(ij),
      ZOC::LARGE => LARGE_ZOC_LUT.ij2i(ij),
    }
  }
  fn ij2j(&self, ij: u64) -> u32 {
    match self {
      ZOC::EMPTY => EMPTY_ZOC.ij2j(ij),
      ZOC::SMALL => SMALL_ZOC_LUT.ij2j(ij),
      ZOC::MEDIUM => MEDIU_ZOC_LUT.ij2j(ij),
      ZOC::LARGE => LARGE_ZOC_LUT.ij2j(ij),
    }
  }
  fn xy2h(&self, x: f64, y: f64) -> u64{
    match self {
      ZOC::EMPTY => EMPTY_ZOC.xy2h(x, y),
      ZOC::SMALL => SMALL_ZOC_LUT.xy2h(x, y),
      ZOC::MEDIUM => MEDIU_ZOC_LUT.xy2h(x, y),
      ZOC::LARGE => LARGE_ZOC_LUT.xy2h(x, y),
    }
  }
}


/// Returns a zorder curve trait implementation according to the given depth.
pub fn get_zoc(depth: u8) -> &'static dyn ZOrderCurve { // on day, It would be possible to add 'const'
  super::super::check_depth(depth);
  #[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
  {
    match depth {
           0  => &EMPTY_ZOC,
       1..=8  => &SMALL_ZOC_BMI, // 2 bits *  8 = 16 bits
       9..=16 => &MEDIU_ZOC_BMI, // 2 bits * 16 = 32 bits
      17..=29 => &LARGE_ZOC_BMI,
      _ => unreachable!(),
    }
  }
  #[cfg(not(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2")))]
    {
      match depth {
             0  => &EMPTY_ZOC,
         1..=8  => &SMALL_ZOC_LUT, // 2 bits *  8 = 16 bits
         9..=16 => &MEDIU_ZOC_LUT, // 2 bits * 16 = 32 bits
        17..=29 => &LARGE_ZOC_LUT,
        _ => unreachable!(),
      }
    }
}


// see config: rustc --print cfg

#[cfg(test)]
mod tests {
  use super::*;
  
  #[test]
  fn test_small(){
    let n = 256u32;
    for i in 0..n {
      for j in 0..n {
        let m1 = SMALL_ZOC_LUT.ij2h(i, j);
        let m2 = SMALL_ZOC_XOR.ij2h(i, j);
        #[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
        let m3 = SMALL_ZOC_BMI.ij2h(i, j);
        assert_eq!(m1, m2);
        #[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
        assert_eq!(m2, m3);
      }
    }
  }

  #[test]
  fn test_mediu(){
    let n = 65536u32;
    for i in (0..n).step_by(7) {
      for j in (0..n).step_by(9) {
        let m1 = MEDIU_ZOC_LUT.ij2h(i, j);
        let m2 = MEDIU_ZOC_XOR.ij2h(i, j);
        #[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
        let m3 = MEDIU_ZOC_BMI.ij2h(i, j);
        // println!("i: {}; j: {}; m1: {}; m2: {}; m3: {}", i, j, m1, m2, m3);
        assert_eq!(m1, m2);
        #[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
        assert_eq!(m2, m3);
      }
    }
  }

  #[test]
  fn test_large(){
    let n = 0xFFFFFFFFu32;
    for i in (0..n).step_by(14999991) {
      for j in (0..n).step_by(1499991) {
        let m1 = LARGE_ZOC_LUT.ij2h(i, j);
        let m2 = LARGE_ZOC_XOR.ij2h(i, j);
        #[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
        let m3 = LARGE_ZOC_BMI.ij2h(i, j);
        assert_eq!(m1, m2);
        #[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
        assert_eq!(m2, m3)
      }
    }
  }
}
