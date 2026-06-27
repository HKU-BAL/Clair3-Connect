use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;

pub fn compute_signal_lengths(record: &Record, is_reverse_strand: bool) -> Option<Vec<i32>> {
    let aux = record.aux(b"mv").ok()?;

    let movements: Vec<i32> = match aux {
        Aux::ArrayI8(a) => a.iter().map(i32::from).collect(),
        Aux::ArrayU8(a) => a.iter().map(i32::from).collect(),
        Aux::ArrayI16(a) => a.iter().map(i32::from).collect(),
        Aux::ArrayU16(a) => a.iter().map(i32::from).collect(),
        Aux::ArrayI32(a) => a.iter().collect(),
        Aux::ArrayU32(a) => a.iter().map(|v| v as i32).collect(),
        _ => return None,
    };
    if movements.len() <= 1 {
        return None;
    }
    let read_len = record.seq_len();
    if read_len == 0 {
        return None;
    }

    let mut signals: Vec<i32> = vec![0; read_len];
    let mut base_index: i64 = -1;

    for &mv in movements.iter().skip(1) {
        if mv != 0 {
            base_index += 1;
            if base_index >= read_len as i64 {
                break;
            }
            signals[base_index as usize] += 1;
        } else {
            if base_index < 0 {
                continue;
            }
            if base_index >= read_len as i64 {
                break;
            }
            signals[base_index as usize] += 1;
        }
    }

    if is_reverse_strand {
        signals.reverse();
    }
    Some(signals)
}

#[cfg(test)]
mod tests {

}
