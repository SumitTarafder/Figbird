import sys

slice= int(sys.argv[3])

with open(sys.argv[2], 'w') as fo:
  with open(sys.argv[1], 'r', encoding="utf8", errors='ignore') as f:
    buffer = ''
    bufferlen_saved=0
    for line in f:
      if line[0] == '>':
          if bufferlen_saved > 0 :
            fo.write(buffer + '\n')
            buffer = ''
          fo.write(line)
          #print(line)
      else:
          buffer+=line[:-1]
      bufferlen=len(buffer)
      bufferlen_saved=bufferlen
      i=0
      while True:
        if bufferlen >= slice:
            fo.write(buffer[slice*i:slice*(i+1)] + '\n')
            bufferlen = bufferlen - slice
            i = i+1
        else:
            fo.write(buffer[slice*i:bufferlen_saved])
            break
      buffer = ''
