// Copyright 2016 Pierre Marijon.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


//! GFF3 format reading and writing.
//!
//! GFF2 definition : http://gmod.org/wiki/GFF2#The_GFF2_File_Format (not yet support)
//! GTF2 definition : http://mblab.wustl.edu/GTF2.html (not yet support)
//! GFF3 definition : http://gmod.org/wiki/GFF3#GFF3_Format
//!
//! # Example
//!
//! ```
//! use std::io;
//! use bio::io::gff;
//! let reader = gff::Reader::new(io::stdin(), gff::GffType::GFF3);
//! ```

use std::io;
use std::fs;
use std::path::Path;
use std::convert::AsRef;
use itertools::Itertools;
use regex::Regex;
use multimap::MultiMap;

use csv;

use utils::{Interval, IntervalError, Strand, StrandError};

/// `GffType`
///
/// We have three format in the GFF family.
/// The change is in the last field of GFF.
/// For each type we have key value separator and field separator
#[derive(Clone, Copy)]
pub enum GffType {
    /// Attribute format is key1=value, key2=value
    GFF3,
    /// Attribute format is key1 value; key2 value
    GFF2,
    /// Same as GFF2 just possible keyword and possible value change
    GTF2,
    /// Any, first field of tuple is key value separator, second is field separator
    Any(u8, u8),
}

impl GffType {
    #[inline]
    fn separator(&self) -> (u8, u8) {
        match *self {
            GffType::GFF3 => (b'=', b','),
            GffType::GFF2 => (b' ', b';'),
            GffType::GTF2 => (b' ', b';'),
            GffType::Any(x, y) => (x, y),
        }
    }
}

/// A GFF reader.
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
    gff_type: GffType,
}

impl Reader<fs::File> {
    /// Read GFF from given file path in given format.
    pub fn from_file<P: AsRef<Path>>(path: P, fileformat: GffType) -> io::Result<Self> {
        fs::File::open(path).map(|f| Reader::new(f, fileformat))
    }
}


impl<R: io::Read> Reader<R> {
    /// Create a new GFF reader given an instance of `io::Read`, in given format.
    pub fn new(reader: R, fileformat: GffType) -> Self {
        Reader {
            inner: csv::Reader::from_reader(reader)
                .delimiter(b'\t')
                .has_headers(false),
            gff_type: fileformat,
        }
    }

    pub fn raw_rows(&mut self) -> RawRows<R> {
        RawRows {
            inner: self.inner.decode(),
        }
    }

    pub fn records(&mut self) -> Records<R> {
        let (delim, term) = self.gff_type.separator();
        let r = format!(
            r" *(?P<key>[^{delim}{term}\t]+){delim}(?P<value>[^{delim}{term}\t]+){term}?",
            delim = delim as char, term = term as char);
        Records {
            inner: self.raw_rows(),
            attributes_re: Regex::new(&r).unwrap(),
        }
    }
}


pub struct RawRows<'a, R: 'a + io::Read> {
    inner: csv::DecodedRecords<'a, R, RawRow>,
}

impl<'a, R: io::Read> Iterator for RawRows<'a, R> {

    type Item = Result<RawRow, GffError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|csvr| csvr.map_err(GffError::from))
    }
}

pub struct Records<'a, R: 'a + io::Read> {
    inner: RawRows<'a, R>,
    attributes_re: Regex,
}

impl<'a, R: io::Read> Iterator for Records<'a, R> {

    type Item = Result<Record, GffError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
            .map(|res| res.and_then(|row| Record::try_from(row, &self.attributes_re)))
    }
}

/// A GFF writer.
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
    gff_type: GffType,
}


impl Writer<fs::File> {
    /// Write to a given file path in given format.
    pub fn to_file<P: AsRef<Path>>(path: P, fileformat: GffType) -> io::Result<Self> {
        fs::File::create(path).map(|f| Writer::new(f, fileformat))
    }
}


impl<W: io::Write> Writer<W> {
    /// Write to a given writer.
    pub fn new(writer: W, fileformat: GffType) -> Self {
        Writer {
            inner: csv::Writer::from_writer(writer)
                .delimiter(b'\t')
                .flexible(true),
            gff_type: fileformat,
        }
    }

    pub fn write(&mut self, record: &Record) -> Result<(), GffError> {
        self.inner
            .encode(record.to_raw_row(&self.gff_type))
            .map_err(GffError::from)
    }
}

quick_error! {

    #[derive(Debug)]
    pub enum GffError {
        ParseFloat(err: ::std::num::ParseFloatError) {
            description(err.description())
            from()
            cause(err)
        }
        ParseInt(err: ::std::num::ParseIntError) {
            description(err.description())
            from()
            cause(err)
        }
        Interval(err: IntervalError) {
            description(err.description())
            from()
            cause(err)
        }
        Strand(err: StrandError) {
            description(err.description())
            from()
            cause(err)
        }
        Csv(err: csv::Error) {
            description(err.description())
            from()
            cause(err)
        }
        InvalidFrame {
            description("GFF record frame must be either 0, 1, 2, if defined.")
        }
    }
}

pub type RawRow = (String, String, String, u64, u64, String, char, char, String);

#[derive(Serialize, Deserialize)]
pub struct Record {
    seqname: String,
    source: Option<String>,
    feature_type: Option<String>,
    interval: Interval<u64>,
    score: Option<f64>,
    strand: Strand,
    frame: Option<u8>,
    attributes: MultiMap<String, String>,
}

pub struct RecordBuilder {
    seqname: String,
    start: u64,
    end: u64,
    source: Option<String>,
    feature_type: Option<String>,
    score: Option<String>,
    strand: Option<char>,
    frame: Option<char>,
    attributes: MultiMap<String, String>,
}

impl RecordBuilder {

    pub fn new<T>(seqname: T, start: u64, end: u64) -> Self
        where T: Into<String>
    {
        RecordBuilder {
            seqname: seqname.into(),
            start: start,
            end: end,
            source: None,
            feature_type: None,
            score: None,
            strand: None,
            frame: None,
            attributes: MultiMap::new(),
        }
    }

    pub fn source<T>(mut self, source: T) -> Self
        where T: Into<String> + AsRef<str>
    {
        self.source = match source.as_ref() {
            "." | "" => None,
            v => Some(v.into()),
        };
        self
    }

    pub fn feature_type<T>(mut self, feature_type: T) -> Self
        where T: Into<String> +  AsRef<str>
    {
        self.feature_type = match feature_type.as_ref() {
            "." | "" => None,
            v => Some(v.into()),
        };
        self
    }

    pub fn score<T>(mut self, raw_score: T) -> Self
        where T: Into<String> + AsRef<str>
    {
        self.score = match raw_score.as_ref() {
            "." | "" => None,
            v => Some(v.into()),
        };
        self
    }

    pub fn strand(mut self, strand_char: char) -> Self {
        self.strand = Some(strand_char);
        self
    }

    pub fn frame(mut self, raw_frame: char) -> Self {
        self.frame = match raw_frame {
            '.' => None,
            v => Some(v),
        };
        self
    }

    pub fn attributes(mut self, attributes: MultiMap<String ,String>) -> Self {
        self.attributes = attributes;
        self
    }

    pub fn build(self) -> Result<Record, GffError> {
        let interval = Interval::new(self.start..self.end)
            .map_err(GffError::from)?;

        let rscore = self.score
            .map(|ref rs| rs.parse::<f64>().map_err(GffError::from));
        let score = match rscore {
            Some(a) => Some(a?),
            None => None,
        };

        let rstrand = self.strand
            .map(|ref v| Strand::from_char(v).map_err(GffError::from));
        let strand = match rstrand {
            Some(a) => a?,
            None => Strand::Unknown,
        };


        let rframe = self.frame
            .map(|rf|
                rf.to_string().parse::<u8>().map_err(GffError::from)
                    .and_then(|num|
                        if num > 2 {
                            Err(GffError::InvalidFrame)
                        } else {
                            Ok(num)
                        }
                    )
            );
        let frame = match rframe {
            Some(a) => Some(a?),
            None => None,
        };

        let record = Record {
            seqname: self.seqname,
            interval: interval,
            source: self.source,
            feature_type: self.feature_type,
            score: score,
            strand: strand,
            frame: frame,
            attributes: self.attributes,
        };
        Ok(record)
    }
}

impl Record {

    fn try_from<'a>(row: RawRow, attributes_re: &Regex) -> Result<Self, GffError> {
        let mut attrs = MultiMap::<String, String>::new();
        let trim_quotes = |s: &str| s.trim_matches('\'').trim_matches('"').to_owned();
        for caps in attributes_re.captures_iter(&row.8) {
            attrs.insert(trim_quotes(&caps["key"]), trim_quotes(&caps["value"]));
        }
        RecordBuilder::new(row.0, row.3 - 1, row.4)
            .source(row.1)
            .feature_type(row.2)
            .score(row.5)
            .strand(row.6)
            .frame(row.7)
            .attributes(attrs)
            .build()
    }

    pub fn to_raw_row(&self, gff_type: &GffType) -> RawRow {
        let (delim, term) = gff_type.separator();
        (self.seqname().to_owned(),
         self.source().unwrap_or(".").to_owned(),
         self.feature_type().unwrap_or(".").to_owned(),
         self.start() + 1,  // start coordinate adjustment
         self.end(),
         self.score().map(|v| v.to_string()).unwrap_or(".".to_owned()),
         match self.strand() {
             Strand::Forward => '+',
             Strand::Reverse => '-',
             Strand::Unknown => '.',
         },
         self.frame().map(|v| (v + 48) as char).unwrap_or('.'),  // make 1u8 into '1', etc.
         self.attributes()
             .iter()
             .map(|(a, b)| format!("{}{}{}", a, delim as char, b))
             .join(::std::str::from_utf8(&[term + 48]).unwrap()),
        )
    }

    /// Sequence name of the feature.
    pub fn seqname(&self) -> &str {
        &self.seqname
    }

    /// Source of the feature.
    pub fn source(&self) -> Option<&str> {
        self.source.as_ref().map(|ref v| v.as_str())
    }

    /// Type of the feature.
    pub fn feature_type(&self) -> Option<&str> {
        self.feature_type.as_ref().map(|ref v| v.as_str())
    }

    /// Start position of feature (0-based).
    pub fn start(&self) -> u64 {
        self.interval.start
    }

    /// End position of feature (0-based, not included).
    pub fn end(&self) -> u64 {
        self.interval.end
    }

    /// Score of feature.
    pub fn score(&self) -> Option<f64> {
        self.score
    }

    /// Strand of the feature.
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Frame of the feature.
    pub fn frame(&self) -> Option<u8> {
        self.frame
    }

    /// Attribute of the feature.
    pub fn attributes(&self) -> &MultiMap<String, String> {
        &self.attributes
    }

    /// Set the seqname of feature.
    pub fn set_seqname<T>(&mut self, seqname: T)
        where T: Into<String>
    {
        self.seqname = seqname.into();
    }

    /// Set the source of feature.
    pub fn set_source<T>(&mut self, source: T)
        where T: Into<String> + AsRef<str>
    {
        self.source = match source.as_ref() {
            "." | "" => None,
            v => Some(v.into()),
        };
    }

    /// Set the type of feature.
    pub fn set_feature_type<T>(&mut self, feature_type: T)
        where T: Into<String> +  AsRef<str>
    {
        self.feature_type = match feature_type.as_ref() {
            "." | "" => None,
            v => Some(v.into()),
        };
    }

    /// Set the start and end coordinate of feature.
    pub fn set_coords(&mut self, start: u64, end: u64) -> Result<(), GffError> {
        let new_interval = Interval::new(start..end).map_err(GffError::from)?;
        self.interval = new_interval;
        Ok(())
    }

    /// Set the score of feature.
    pub fn set_score(&mut self, score: f64) {
        self.score = Some(score);
    }

    /// Set the strand of feature.
    pub fn set_strand(&mut self, strand: Strand) {
        self.strand = strand;
    }

    /// Set the frame of feature.
    pub fn set_frame(&mut self, frame: u8) -> Result<(), GffError> {
        if frame > 2 {
            Err(GffError::InvalidFrame)
        } else {
            self.frame = Some(frame);
            Ok(())
        }
    }

    /// Set the attributes of feature.
    pub fn set_attributes(&mut self, attributes: MultiMap<String, String>) {
        self.attributes = attributes;
    }

    /// Get mutable reference on attributes of feature.
    pub fn attributes_mut(&mut self) -> &mut MultiMap<String, String> {
        &mut self.attributes
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use utils::Strand;
    use multimap::MultiMap;

    const GFF_FILE: &'static [u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote=Removed,ID=test
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tNote=ATP-dependent protease subunit HslV,ID=PRO_0000148105
";
    //required because MultiMap iter on element randomly
    const GFF_FILE_ONE_ATTRIB: &'static [u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote=Removed
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tID=PRO_0000148105
";

    const GTF_FILE: &'static [u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote Removed;ID test
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tNote ATP-dependent;ID PRO_0000148105
";

    // Another variant of GTF file, modified from a published GENCODE GTF file.
    const GTF_FILE_2: &'static [u8] = b"chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; transcript_id \"ENST00000456328.2\"; gene_type \"transcribed_unprocessed_pseudogene\"";

    // GTF file with duplicate attribute keys, taken from a published GENCODE GTF file.
    const GTF_FILE_DUP_ATTR_KEYS: &'static [u8] = b"chr1\tENSEMBL\ttranscript\t182393\t184158\t.\t+\t.\tgene_id \"ENSG00000279928.1\"; transcript_id \"ENST00000624431.1\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"FO538757.2\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"FO538757.2-201\"; level 3; protein_id \"ENSP00000485457.1\"; transcript_support_level \"1\"; tag \"basic\"; tag \"appris_principal_1\";";

    //required because MultiMap iter on element randomly
    const GTF_FILE_ONE_ATTRIB: &'static [u8] = b"P0A7B8\tUniProtKB\tInitiator methionine\t1\t1\t.\t.\t.\tNote Removed
P0A7B8\tUniProtKB\tChain\t2\t176\t50\t+\t.\tID PRO_0000148105
";

    #[test]
    fn test_raw_rows() {
        let seqname = ["P0A7B8", "P0A7B8"];
        let source = ["UniProtKB", "UniProtKB"];
        let feature_type = ["Initiator methionine", "Chain"];
        let starts = [1, 2];
        let ends = [1, 176];
        let scores = [".", "50"];
        let strand = ['.', '+'];
        let frame = ['.', '.'];
        let attributes = ["Note=Removed,ID=test",
                          "Note=ATP-dependent protease subunit HslV,ID=PRO_0000148105"];

        let mut reader = Reader::new(GFF_FILE, GffType::GFF3);
        for (i, r) in reader.raw_rows().enumerate() {
            let row = r.unwrap();
            assert_eq!(row.0, seqname[i]);
            assert_eq!(row.1, source[i]);
            assert_eq!(row.2, feature_type[i]);
            assert_eq!(row.3, starts[i]);
            assert_eq!(row.4, ends[i]);
            assert_eq!(row.5, scores[i]);
            assert_eq!(row.6, strand[i]);
            assert_eq!(row.7, frame[i]);
            assert_eq!(row.8, attributes[i]);
        }
    }

    #[test]
    fn test_reader_gff3() {
        let seqname = ["P0A7B8", "P0A7B8"];
        let source = [Some("UniProtKB"), Some("UniProtKB")];
        let feature_type = [Some("Initiator methionine"), Some("Chain")];
        let starts = [0, 1];
        let ends = [1, 176];
        let scores = [None, Some(50f64)];
        let strand = [Strand::Unknown, Strand::Forward];
        let frame = [None, None];
        let mut attributes = [MultiMap::new(), MultiMap::new()];
        attributes[0].insert("ID".to_owned(), "test".to_owned());
        attributes[0].insert("Note".to_owned(), "Removed".to_owned());
        attributes[1].insert("ID".to_owned(), "PRO_0000148105".to_owned());
        attributes[1].insert("Note".to_owned(),
                             "ATP-dependent protease subunit HslV".to_owned());

        let mut reader = Reader::new(GFF_FILE, GffType::GFF3);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
            assert_eq!(record.seqname(), seqname[i]);
            assert_eq!(record.source(), source[i]);
            assert_eq!(record.feature_type(), feature_type[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.score(), scores[i]);
            assert!(record.strand().is_unknown() || record.strand() == strand[i]);
            assert_eq!(record.frame(), frame[i]);
            assert_eq!(record.attributes(), &attributes[i]);
        }
    }

    #[test]
    fn test_reader_gtf2() {
        let seqname = ["P0A7B8", "P0A7B8"];
        let source = [Some("UniProtKB"), Some("UniProtKB")];
        let feature_type = [Some("Initiator methionine"), Some("Chain")];
        let starts = [0, 1];
        let ends = [1, 176];
        let scores = [None, Some(50f64)];
        let strand = [Strand::Unknown, Strand::Forward];
        let frame = [None, None];
        let mut attributes = [MultiMap::new(), MultiMap::new()];
        attributes[0].insert("ID".to_owned(), "test".to_owned());
        attributes[0].insert("Note".to_owned(), "Removed".to_owned());
        attributes[1].insert("ID".to_owned(), "PRO_0000148105".to_owned());
        attributes[1].insert("Note".to_owned(), "ATP-dependent".to_owned());

        let mut reader = Reader::new(GTF_FILE, GffType::GTF2);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
            assert_eq!(record.seqname(), seqname[i]);
            assert_eq!(record.source(), source[i]);
            assert_eq!(record.feature_type(), feature_type[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.score(), scores[i]);
            assert!(record.strand().is_unknown() || record.strand() == strand[i]);
            assert_eq!(record.frame(), frame[i]);
            assert_eq!(record.attributes(), &attributes[i]);
        }
    }

    #[test]
    fn test_reader_gtf2_2() {
        let seqname = ["chr1", "chr1"];
        let source = [Some("HAVANA"), Some("HAVANA")];
        let feature_type = [Some("gene"), Some("transcript")];
        let starts = [11868, 11868];
        let ends = [14409, 14409];
        let scores = [None, None];
        let strand = [Strand::Forward, Strand::Forward];
        let frame = [None, None];
        let mut attributes = [MultiMap::new(), MultiMap::new()];
        attributes[0].insert("gene_id".to_owned(), "ENSG00000223972.5".to_owned());
        attributes[0].insert("gene_type".to_owned(),
                             "transcribed_unprocessed_pseudogene".to_owned());
        attributes[1].insert("gene_id".to_owned(), "ENSG00000223972.5".to_owned());
        attributes[1].insert("transcript_id".to_owned(), "ENST00000456328.2".to_owned());
        attributes[1].insert("gene_type".to_owned(),
                             "transcribed_unprocessed_pseudogene".to_owned());

        let mut reader = Reader::new(GTF_FILE_2, GffType::GTF2);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();
            assert_eq!(record.seqname(), seqname[i]);
            assert_eq!(record.source(), source[i]);
            assert_eq!(record.feature_type(), feature_type[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.score(), scores[i]);
            assert!(record.strand().is_unknown() || record.strand() == strand[i]);
            assert_eq!(record.frame(), frame[i]);
            assert_eq!(record.attributes(), &attributes[i]);
        }
    }

    #[test]
    fn test_reader_gtf2_dup_attr_keys() {
        let mut reader = Reader::new(GTF_FILE_DUP_ATTR_KEYS, GffType::GTF2);
        let mut records = reader.records().collect::<Vec<_>>();
        assert_eq!(records.len(), 1);
        let record = records.pop().unwrap().expect("expected one record");
        assert_eq!(record.attributes.get("tag"),
                   Some(&"basic".to_owned()));
        assert_eq!(record.attributes.get_vec("tag"),
                   Some(&vec!["basic".to_owned(), "appris_principal_1".to_owned()]));
    }

    #[test]
    fn test_writer_gff3() {
        let mut reader = Reader::new(GFF_FILE_ONE_ATTRIB, GffType::GFF3);
        let mut writer = Writer::new(vec![], GffType::GFF3);
        for r in reader.records() {
            writer
                .write(&r.ok().expect("Error reading record"))
                .ok()
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.as_string(),
                   String::from_utf8_lossy(GFF_FILE_ONE_ATTRIB))
    }

    #[test]
    fn test_writer_gtf2() {
        let mut reader = Reader::new(GTF_FILE_ONE_ATTRIB, GffType::GTF2);
        let mut writer = Writer::new(vec![], GffType::GTF2);
        for r in reader.records() {
            writer
                .write(&r.ok().expect("Error reading record"))
                .ok()
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.as_string(),
                   String::from_utf8_lossy(GTF_FILE_ONE_ATTRIB))
    }

    #[test]
    fn test_convert_gtf2_to_gff3() {
        let mut reader = Reader::new(GTF_FILE_ONE_ATTRIB, GffType::GTF2);
        let mut writer = Writer::new(vec![], GffType::GFF3);
        for r in reader.records() {
            writer
                .write(&r.ok().expect("Error reading record"))
                .ok()
                .expect("Error writing record");
        }
        assert_eq!(writer.inner.as_string(),
                   String::from_utf8_lossy(GFF_FILE_ONE_ATTRIB))
    }
}
