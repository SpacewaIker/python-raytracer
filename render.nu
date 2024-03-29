
def main [infile: string, num_tasks: int, outfile: string] {
  pueued --daemonize
  pueue parallel $num_tasks

  mkdir temp

  mut tids = []

  for n in ..<$num_tasks {
    let tid = (pueue add -p -- python provided/main.py --infile $infile --outfile $"temp/($n)" --subimage $n --tasks $num_tasks)
    $tids ++= $tid
  }

  let tid = (pueue add -p -a ...$tids -- python provided/glue.py --dir temp --outfile $outfile)

  pueue follow $tid

  rm -r temp

  pueue clean
  pueue shutdown
}

