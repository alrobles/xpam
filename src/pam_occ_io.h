// pam_occ_io.h -- occurrence-pair sources for the PAM engine.
//
// The compute loop consumes (cell, species) pairs sorted by cell exactly
// once, in a single merge pass against raster strips.  OccSource abstracts
// where those pairs live:
//
//   MemOccSource  -- vectors already in memory (index array sorted by cell)
//   FileOccSource -- chunked stream from a sorted binary file, so the pair
//                    table never materializes in R *or* C++ (Phase C)
//
// File format ("XPOC" v1, columnar, written by xpam_write_occ_bin() in R):
//   offset  0:  char[4]  "XPOC"
//   offset  4:  uint32   version = 1
//   offset  8:  uint64   n_pairs
//   offset 16:  uint64   n_species
//   offset 24:  float64  cells[n_pairs]   (0-based, sorted ascending)
//   offset 24 + 8*n: int32 sps[n_pairs]   (0-based species codes)
//
// cells are stored as float64 because R cannot emit raw int64; every raster
// cell index is < 2^53 so the representation is exact.
//
// Chunk size bounds the streaming window (~48 MiB at the default 4M pairs).

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace pam {

class OccSource {
public:
  virtual ~OccSource() = default;
  virtual std::int64_t n_pairs()   const = 0;
  virtual std::int64_t n_species() const = 0;
  // Yield the next (cell, sp) in cell-sorted order; false when exhausted.
  virtual bool next(std::int64_t& cell, std::int32_t& sp) = 0;
  // Inputs skipped as unusable (not out-of-range pairs): default none.
  virtual std::int64_t n_skipped() const { return 0; }
};

// In-memory source over index arrays, sorted by cell on first use.
class MemOccSource : public OccSource {
public:
  MemOccSource(const int* sp, const int* cell, std::int64_t n,
               std::int64_t n_species, std::int64_t ncells,
               std::int64_t* n_oor)
      : sp_(sp), cell_(cell), n_(n), nsp_(n_species) {
    order_.reserve((std::size_t)n_);
    for (std::int64_t p = 0; p < n_; ++p) {
      const int c = cell_[p], s = sp_[p];
      if (c == NA_INT32 || s == NA_INT32) continue;
      if (c < 0 || c >= ncells || s < 0 || s >= nsp_) { ++(*n_oor); continue; }
      order_.push_back(((std::uint64_t)(std::uint32_t)c << 32)
                       | (std::uint32_t)p);
    }
    std::sort(order_.begin(), order_.end());
  }
  std::int64_t n_pairs()   const override { return (std::int64_t)order_.size(); }
  std::int64_t n_species() const override { return nsp_; }
  bool next(std::int64_t& cell, std::int32_t& sp) override {
    if (i_ >= order_.size()) return false;
    const std::uint64_t key = order_[i_++];
    cell = (std::int64_t)(key >> 32);
    sp   = sp_[key & 0xffffffffu];
    return true;
  }
private:
  const int *sp_, *cell_;
  std::int64_t n_, nsp_;
  std::vector<std::uint64_t> order_;
  std::size_t i_ = 0;
  static constexpr int NA_INT32 = std::numeric_limits<int>::min();
};

// k-way merge of cell-sorted sources (one per polygon dataset, e.g. the
// per-group .shp files of an IUCN-style collection).  A min-heap on the
// head pair of each source preserves the global (cell, sp) order the
// engine's merge pass requires; equal pairs emitted by overlapping
// datasets are deduplicated on the fly.  Species codes must already be
// consistent across sources — the caller passes one shared level vector
// to every child.
class MergeOccSource : public OccSource {
public:
  explicit MergeOccSource(std::vector<std::unique_ptr<OccSource>> srcs)
      : srcs_(std::move(srcs)) {
    heap_.reserve(srcs_.size());
    for (std::size_t k = 0; k < srcs_.size(); ++k) {
      std::int64_t c; std::int32_t s;
      if (srcs_[k]->next(c, s))
        heap_.push_back(Node{c, s, k});
      np_ += srcs_[k]->n_pairs() < 0 ? 0 : srcs_[k]->n_pairs();
      any_unknown_ |= srcs_[k]->n_pairs() < 0;
    }
    std::make_heap(heap_.begin(), heap_.end(), Greater{});
  }
  MergeOccSource(const MergeOccSource&) = delete;
  MergeOccSource& operator=(const MergeOccSource&) = delete;

  // Sum of children, or -1 while any child cannot know (polygon sources).
  std::int64_t n_pairs()   const override { return any_unknown_ ? -1 : np_; }
  std::int64_t n_species() const override {
    return srcs_.empty() ? 0 : srcs_[0]->n_species();
  }
  std::int64_t n_skipped() const override {
    std::int64_t n = 0;
    for (const auto& s : srcs_) n += s->n_skipped();
    return n;
  }

  bool next(std::int64_t& cell, std::int32_t& sp) override {
    while (!heap_.empty()) {
      std::pop_heap(heap_.begin(), heap_.end(), Greater{});
      const Node n = heap_.back(); heap_.pop_back();
      std::int64_t c2; std::int32_t s2;
      if (srcs_[n.k]->next(c2, s2)) {          // refill the popped source
        heap_.push_back(Node{c2, s2, n.k});
        std::push_heap(heap_.begin(), heap_.end(), Greater{});
      }
      if (n.cell == last_cell_ && n.sp == last_sp_) continue;  // cross-source dup
      last_cell_ = n.cell; last_sp_ = n.sp;
      cell = n.cell; sp = n.sp;
      return true;
    }
    return false;
  }

private:
  struct Node { std::int64_t cell; std::int32_t sp; std::size_t k; };
  struct Greater {
    bool operator()(const Node& a, const Node& b) const {
      return a.cell != b.cell ? a.cell > b.cell : a.sp > b.sp;
    }
  };
  std::vector<std::unique_ptr<OccSource>> srcs_;
  std::vector<Node> heap_;
  std::int64_t np_ = 0, last_cell_ = -1;
  std::int32_t last_sp_ = -1;
  bool any_unknown_ = false;
};

// Chunked file source; the file must already be cell-sorted (writer's
// contract).  Sortedness is verified incrementally while draining.
class FileOccSource : public OccSource {
public:
  FileOccSource(const std::string& path, std::size_t chunk_pairs = 1u << 22,
                std::int64_t* n_oor = nullptr, std::int64_t ncells = -1) {
    f_ = std::fopen(path.c_str(), "rb");
    if (!f_) throw std::runtime_error("cannot open occurrence file: " + path);
    char magic[4];
    std::uint32_t ver;
    std::uint64_t n64, nsp64;
    if (std::fread(magic, 1, 4, f_) != 4 || std::memcmp(magic, "XPOC", 4) ||
        std::fread(&ver, 4, 1, f_) != 1 || ver != 1 ||
        std::fread(&n64, 8, 1, f_) != 1 ||
        std::fread(&nsp64, 8, 1, f_) != 1) {
      std::fclose(f_);
      throw std::runtime_error("bad XPOC header in: " + path);
    }
    n_ = (std::int64_t)n64; nsp_ = (std::int64_t)nsp64;
    cells_off_ = std::ftell(f_);
    sps_off_   = cells_off_ + n_ * 8;
    chunk_ = chunk_pairs ? chunk_pairs : 1;
    ncells_ = ncells; oor_ = n_oor;
  }
  ~FileOccSource() override { if (f_) std::fclose(f_); }
  FileOccSource(const FileOccSource&) = delete;
  FileOccSource& operator=(const FileOccSource&) = delete;

  std::int64_t n_pairs()   const override { return n_; }
  std::int64_t n_species() const override { return nsp_; }

  bool next(std::int64_t& cell, std::int32_t& sp) override {
    if (pos_ >= n_) return false;
    if (i_ >= end_) refill();
    cell = cells_[i_];
    sp   = sps_[i_++];
    if (cell < last_cell_) {
      std::fclose(f_); f_ = nullptr;
      throw std::runtime_error("occurrence file is not cell-sorted");
    }
    last_cell_ = cell;
    ++pos_;
    if (ncells_ >= 0 && (cell < 0 || cell >= ncells_ || sp < 0 || sp >= nsp_)) {
      if (oor_) ++(*oor_);
      return next(cell, sp);           // skip out-of-range, keep streaming
    }
    return true;
  }

private:
  void refill() {
    const std::size_t k = (std::size_t)std::min<std::int64_t>(
        (std::int64_t)chunk_, n_ - pos_);
    cells_.resize(k); sps_.resize(k);
    if (std::fseek(f_, cells_off_ + pos_ * 8, SEEK_SET) ||
        std::fread(cells_.data(), 8, k, f_) != k ||
        std::fseek(f_, sps_off_ + pos_ * 4, SEEK_SET) ||
        std::fread(sps_.data(), 4, k, f_) != k) {
      std::fclose(f_); f_ = nullptr;
      throw std::runtime_error("truncated occurrence file");
    }
    i_ = 0; end_ = k;
  }

  std::FILE* f_ = nullptr;
  std::int64_t n_ = 0, nsp_ = 0, pos_ = 0, ncells_ = -1;
  std::int64_t cells_off_ = 0, sps_off_ = 0;
  std::size_t chunk_ = 0, i_ = 0, end_ = 0;
  std::int64_t last_cell_ = -1;
  std::vector<double>   cells_;
  std::vector<std::int32_t> sps_;
  std::int64_t* oor_ = nullptr;
};

// Streaming writer for the XPOC v1 layout.  The header stores n_pairs,
// which is only known once the stream ends, so pairs go to two temp files
// (cells f64 / sps i32, appended in emission order -- the caller's cell
// sort) and finish() splices them into the columnar file.
class XpocWriter {
public:
  XpocWriter(const std::string& path, std::int64_t n_species,
             std::size_t buf_pairs = 1u << 20)
      : path_(path), nsp_(n_species) {
    fc_ = std::fopen((path + ".cells.tmp").c_str(), "wb");
    fs_ = std::fopen((path + ".sps.tmp").c_str(), "wb");
    if (!fc_ || !fs_) {
      if (fc_) std::fclose(fc_);
      if (fs_) std::fclose(fs_);
      throw std::runtime_error("cannot open XPOC temp files for: " + path);
    }
    cbuf_.reserve(buf_pairs); sbuf_.reserve(buf_pairs);
  }
  ~XpocWriter() {
    try { finish(); } catch (...) {}   // best effort; explicit finish() throws
  }
  XpocWriter(const XpocWriter&) = delete;
  XpocWriter& operator=(const XpocWriter&) = delete;

  void add(std::int64_t cell, std::int32_t sp) {
    cbuf_.push_back((double)cell); sbuf_.push_back(sp);
    if (cbuf_.size() == cbuf_.capacity()) flush_();
    ++n_;
  }
  std::int64_t n_pairs() const { return n_; }

  // Drop the temp files without producing the final artefact (used when
  // the stream was interrupted); a no-op after finish().
  void abandon() {
    if (done_) return;
    done_ = true;
    if (fc_) { std::fclose(fc_); fc_ = nullptr; }
    if (fs_) { std::fclose(fs_); fs_ = nullptr; }
    std::remove((path_ + ".cells.tmp").c_str());
    std::remove((path_ + ".sps.tmp").c_str());
  }

  void finish() {
    if (done_) return;
    done_ = true;
    flush_();
    std::fclose(fc_); std::fclose(fs_); fc_ = fs_ = nullptr;
    std::FILE* out = std::fopen(path_.c_str(), "wb");
    if (!out)
      throw std::runtime_error("cannot open XPOC output: " + path_);
    std::fwrite("XPOC", 1, 4, out);
    const std::uint32_t ver = 1;
    const std::uint64_t n64 = (std::uint64_t)n_, nsp64 = (std::uint64_t)nsp_;
    std::fwrite(&ver, 4, 1, out);
    std::fwrite(&n64, 8, 1, out);
    std::fwrite(&nsp64, 8, 1, out);
    copy_(path_ + ".cells.tmp", out, 8);
    copy_(path_ + ".sps.tmp", out, 4);
    if (std::fclose(out) != 0)
      throw std::runtime_error("failed writing XPOC output: " + path_);
    std::remove((path_ + ".cells.tmp").c_str());
    std::remove((path_ + ".sps.tmp").c_str());
  }

private:
  void flush_() {
    if (cbuf_.empty()) return;
    std::fwrite(cbuf_.data(), 8, cbuf_.size(), fc_);
    std::fwrite(sbuf_.data(), 4, sbuf_.size(), fs_);
    cbuf_.clear(); sbuf_.clear();
  }
  static void copy_(const std::string& from, std::FILE* to,
                    std::size_t elem) {
    std::FILE* f = std::fopen(from.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot reopen XPOC temp: " + from);
    std::vector<char> buf(1u << 20);
    std::size_t got;
    while ((got = std::fread(buf.data(), elem, buf.size() / elem, f)) > 0)
      std::fwrite(buf.data(), elem, got, to);
    std::fclose(f);
  }

  std::string path_;
  std::FILE *fc_ = nullptr, *fs_ = nullptr;
  std::vector<double>   cbuf_;
  std::vector<std::int32_t> sbuf_;
  std::int64_t n_ = 0, nsp_ = 0;
  bool done_ = false;
};

// Tees the emitted pair stream into an XPOC file -- the on-disk sparse
// PAM.  Applied as a decorator *after* source construction (including a
// MergeOccSource), so what lands on disk is exactly what the engine
// consumes: merged and deduplicated, not a per-dataset fragment.
class DumpOccSource : public OccSource {
public:
  DumpOccSource(std::unique_ptr<OccSource> inner, const std::string& path)
      : inner_(std::move(inner)) {
    dump_.reset(new XpocWriter(path, inner_->n_species()));
  }
  ~DumpOccSource() override {
    if (dump_) dump_->abandon();    // no-op once finish() ran
  }
  DumpOccSource(const DumpOccSource&) = delete;
  DumpOccSource& operator=(const DumpOccSource&) = delete;

  std::int64_t n_pairs()   const override { return inner_->n_pairs(); }
  std::int64_t n_species() const override { return inner_->n_species(); }
  std::int64_t n_skipped() const override { return inner_->n_skipped(); }

  bool next(std::int64_t& cell, std::int32_t& sp) override {
    if (!inner_->next(cell, sp)) {
      if (dump_) { dump_->finish(); dump_.reset(); }
      return false;
    }
    dump_->add(cell, sp);
    return true;
  }

private:
  std::unique_ptr<OccSource> inner_;
  std::unique_ptr<XpocWriter> dump_;
};

} // namespace pam
