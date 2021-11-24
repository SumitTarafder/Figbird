import sys
with open(sys.argv[2], 'w') as fo:
  with open(sys.argv[1], 'r') as f:
    buffer = ''
    for line in f:
      if line[0] == '>':
          if len(buffer) > 0 :
            fo.write(buffer + '\n')
            buffer = ''
          fo.write(line)
      else:
          buffer+=line[:-1]
      while len(buffer) >= 60:
            fo.write(buffer[:60] + '\n')
            buffer = buffer[60:]

    fo.write(buffer)
