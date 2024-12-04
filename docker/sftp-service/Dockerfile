FROM atmoz/sftp:latest

# Set environment variables for SFTP server
ENV SFTP_USERS=user:password:1001

# Expose the SFTP port
EXPOSE 22

# Create a directory for SFTP files
RUN mkdir -p /home/user/uploads