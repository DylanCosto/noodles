//! CRAM container blocks.

mod compression_method;
mod content_type;

pub use self::{compression_method::CompressionMethod, content_type::ContentType};

/// A block content ID.
pub type ContentId = i32;
